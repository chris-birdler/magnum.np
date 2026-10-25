#!/usr/bin/env python
"""
GPU profile of the time step of a production run (PROTOKOLL 3.36, audit items A1 ... A5): does the GPU (kernel time)
or the CPU side (Python, kernel launches) limit a stage? The wrapper runs run_loops.main with the given arguments (the
code path of production). One call of LLGSolver.step integrates one sample interval. After --skip sample intervals
(warm-up, torch.compile) it measures --plain sample intervals without the profiler (wall clock; the device is
synchronized only at the two ends), then --window sample intervals with torch.profiler (CPU and CUDA activities).
In the profiled window each field term h(), each stage (LLGSolver.dm) and each LLGSolver.step carries a
record_function label: the GPU time of the kernels is split into the terms, the rest of the stage (field sum, LLG
right-hand side) and the rest of the step (RKF45 combination, error norm, step control). Then the run stops.

    CUDA_DEVICE=0 python prof_step.py --skip 300 --plain 40 --window 10 --json out.json -- <run_loops arguments>
    CUDA_DEVICE=0 python prof_step.py --skip 30 --plain 96 --window 0 --kcache DIR --json out.json -- <arguments>

Speed measurement of PLAN 14 (--window 0, no profiler): ms per sample interval, per RKF45 stage and tries per sample
interval of the plain window and of its two halves, the set-up time (start of this script to the first sample
interval, with the kernel build if the cache --kcache does not have the kernel) and the peak GPU memory.
--trace records after each sample interval t, the count of RKF45 tries, <Ms m_x> and the count of accepted steps
(state._step), for check G of PLAN 14a/14b. --null reverse|rotate makes a null run as step_dump.py --k1_null (the
field terms summed in the reverse or in the rotated order): the scale of the rounding for G.

GPU busy share = GPU time per stage (profiled window) / wall time per stage (plain window). Near 1: the GPU limits the
stage. Well below 1: the CPU side limits it.
"""
import argparse
import json
import pathlib
import sys
import time

import torch
from torch.autograd import DeviceType
from torch.profiler import ProfilerActivity, profile, record_function

ST = {"steps": 0, "dm": 0, "tries": 0, "phase": "warm", "wrapped": set()}
T0 = time.perf_counter()
A = None


def sync():
    if torch.cuda.is_available():
        torch.cuda.synchronize()


def labelled(cls, meth, name, count=None):
    """Wrap cls.meth: count its calls (key count) and, in the profiled window, give it a record_function label."""
    if (cls, meth) in ST["wrapped"]:
        return
    ST["wrapped"].add((cls, meth))
    f = getattr(cls, meth)

    def g(self, *a, **k):
        if count is not None:
            ST[count] += 1
        if ST["phase"] == "prof":
            with record_function(name):
                return f(self, *a, **k)
        return f(self, *a, **k)

    setattr(cls, meth, g)


def analyse(prof, wall, n_dm, n_try, n_steps):
    ev = prof.events()
    cpu = [e for e in ev if e.device_type == DeviceType.CPU]
    # device events: kernels, copies and fills; a record_function range (also the ranges of torch.compile) has a
    # GPU-side event too, marked as a user annotation, and it is not GPU work
    dev = [e for e in ev if e.device_type == DeviceType.CUDA and not e.is_user_annotation]
    iv = sorted((e.time_range.start, e.time_range.end) for e in dev)
    busy, cur_s, cur_e = 0.0, None, None                      # union of the kernel intervals (us)
    for s, e in iv:
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                busy += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        busy += cur_e - cur_s
    ksum = sum(e.time_range.elapsed_us() for e in dev)
    # GPU time of each label without its inner labels. A kernel (or copy) has the correlation id of its launch call
    # on the CPU (cudaLaunchKernel, cuLaunchKernel for the Triton kernels of torch.compile, cudaMemcpyAsync ...);
    # the kernel goes to the innermost label on the call chain of that launch call
    labels = {}
    for e in cpu:
        if e.name.startswith("L:"):
            r = labels.setdefault(e.name[2:], {"calls": 0, "gpu_ms_excl": 0.0, "cpu_ms_incl": 0.0})
            r["calls"] += 1
            r["cpu_ms_incl"] += e.time_range.elapsed_us() / 1e3
    launch = {e.id: e for e in cpu if e.name.startswith("cu")}
    kern = {}
    for d in dev:
        p = launch.get(d.id)
        while p is not None and not p.name.startswith("L:"):
            p = p.cpu_parent
        lab = p.name[2:] if p is not None else "(no label)"
        r = labels.setdefault(lab, {"calls": 0, "gpu_ms_excl": 0.0, "cpu_ms_incl": 0.0})
        r["gpu_ms_excl"] += d.time_range.elapsed_us() / 1e3
        k = kern.setdefault(d.name, [0, 0.0, {}])
        k[0] += 1
        k[1] += d.time_range.elapsed_us()
        k[2][lab] = k[2].get(lab, 0.0) + d.time_range.elapsed_us() / 1e3
    cpu_ops = sorted(((e.key, e.count, e.self_cpu_time_total / 1e3) for e in prof.key_averages()
                      if not e.key.startswith("L:")), key=lambda r: -r[2])[:30]
    return {"wall_s": wall, "stages": n_dm, "tries": n_try, "sample_intervals": n_steps,
            "gpu_busy_ms": busy / 1e3, "gpu_kernel_sum_ms": ksum / 1e3, "kernels": len(dev),
            "labels": labels,
            "top_kernels": sorted(([k, v[0], v[1] / 1e3, v[2]] for k, v in kern.items()), key=lambda r: -r[2]),
            "top_cpu_self": [list(r) for r in cpu_ops]}


def report_plain(res):
    p = res["plain"]
    for k, v in [("all", p)] + list(zip(("first half", "second half"), p["halves"])):
        v["ms_per_sample"] = 1e3 * v["wall_s"] / v["sample_intervals"]
        v["ms_per_stage"] = 1e3 * v["wall_s"] / v["stages"]
        v["tries_per_sample"] = v["tries"] / v["sample_intervals"]
        print("[prof_step] plain %-11s: %3d sample intervals, %6d stages, %5d tries, %8.2f s: %8.2f ms per sample "
              "interval, %.4f ms per stage, %.2f tries per sample interval"
              % (k, v["sample_intervals"], v["stages"], v["tries"], v["wall_s"], v["ms_per_sample"],
                 v["ms_per_stage"], v["tries_per_sample"]), flush=True)
    h = p["halves"]
    res["halves_rel_diff"] = abs(h[1]["ms_per_stage"] / h[0]["ms_per_stage"] - 1.0)
    print("[prof_step] set-up %.1f s; peak GPU memory %.2f GiB allocated, %.2f GiB reserved; ms per stage of the two "
          "halves differ by %.1f %%" % (res["setup_s"], res["peak_alloc_GiB"], res["peak_reserved_GiB"],
                                        100 * res["halves_rel_diff"]), flush=True)
    if A.json:
        pathlib.Path(A.json).write_text(json.dumps(res, indent=1))


def report(res):
    p, w = res["plain"], res["window"]
    ms_stage = 1e3 * p["wall_s"] / p["stages"]
    gpu_stage = w["gpu_busy_ms"] / w["stages"]
    res["ms_per_stage_plain"] = ms_stage
    res["gpu_ms_per_stage"] = gpu_stage
    res["gpu_busy_share"] = gpu_stage / ms_stage
    res["kernels_per_stage"] = w["kernels"] / w["stages"]
    res["stages_per_sample"] = p["stages"] / p["sample_intervals"]
    res["stages_per_try"] = p["stages"] / max(p["tries"], 1)
    print("[prof_step] plain: %d sample intervals, %d stages, %d RKF45 tries, %.2f s wall: %.3f ms per stage"
          % (p["sample_intervals"], p["stages"], p["tries"], p["wall_s"], ms_stage), flush=True)
    print("[prof_step] profiled: %d stages, %.2f s wall (with profiler), GPU busy %.1f ms (%.3f ms per stage), "
          "%d kernels (%.1f per stage)" % (w["stages"], w["wall_s"], w["gpu_busy_ms"], gpu_stage, w["kernels"],
                                           res["kernels_per_stage"]), flush=True)
    print("[prof_step] GPU busy share of the plain wall time per stage: %.3f" % res["gpu_busy_share"], flush=True)
    print("[prof_step] GPU time attributed to labels: %.1f ms of %.1f ms kernel time"
          % (sum(v["gpu_ms_excl"] for v in w["labels"].values()), w["gpu_kernel_sum_ms"]), flush=True)
    print("[prof_step] %-34s %8s %22s %26s" % ("label", "calls", "GPU ms/stage (own)", "CPU ms/stage (with inner)"))
    for k, v in sorted(w["labels"].items(), key=lambda kv: -kv[1]["gpu_ms_excl"]):
        print("[prof_step] %-34s %8d %22.4f %26.4f" % (k, v["calls"], v["gpu_ms_excl"] / w["stages"],
                                                       v["cpu_ms_incl"] / w["stages"]))
    print("[prof_step] top kernels (name, calls, GPU ms per stage):")
    for k, n, t, lab in w["top_kernels"][:20]:
        print("[prof_step]   %-80s %7d %9.4f  %s" % (k[:80], n, t / w["stages"], max(lab, key=lab.get)))
    if A.json:
        pathlib.Path(A.json).write_text(json.dumps(res, indent=1))


def main(argv):
    global A
    if "--" not in argv:
        raise SystemExit("usage: prof_step.py [options] -- <run_loops arguments>")
    i = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip", type=int, default=300, help="sample intervals before the measurement (warm-up)")
    ap.add_argument("--plain", type=int, default=40, help="sample intervals of the wall-clock window")
    ap.add_argument("--window", type=int, default=10, help="sample intervals of the profiled window")
    ap.add_argument("--json", default=None)
    ap.add_argument("--kcache", default=None, help="kernel cache folder (kcache.py; one for each geometry and precision)")
    ap.add_argument("--trace", action="store_true",
                    help="record t, the count of RKF45 tries and <Ms m_x> after each sample interval (check G of "
                         "PLAN 14a; one reduction for each sample interval)")
    ap.add_argument("--null", choices=["reverse", "rotate"], default=None,
                    help="null run as step_dump.py --k1_null (scale of the rounding for check G; not the kernel null "
                         "run: 1 + 2^-52 rounds to 1 in single precision)")
    ap.add_argument("--rl_dir", default=".", help="folder of run_loops.py (default: the current folder)")
    A = ap.parse_args(argv[:i])
    rl = argv[i + 1:]
    sys.path.insert(0, str(pathlib.Path(A.rl_dir).resolve()))
    import magnumnp
    from magnumnp.solvers.ode_solvers.rkf45 import RKF45
    import run_loops
    if A.kcache:
        sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
        import kcache
        kcache.use(A.kcache)
    if A.null:
        sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
        from step_dump import k1_null
        k1_null(A.null)
        print("[prof_step] null run: %s" % A.null, flush=True)
    res = {"argv": rl, "skip": A.skip, "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu",
           "torch": torch.__version__, "null": A.null}
    labelled(magnumnp.LLGSolver, "dm", "L:stage (LLGSolver.dm)", count="dm")
    labelled(RKF45, "_try_step", "L:RKF45 try", count="tries")
    f_step = magnumnp.LLGSolver.step
    c = {}

    def step(self, state, dt, **k):
        if ST["steps"] == 0:
            for term in self._terms:                            # the field terms of this run
                labelled(type(term), "h", "L:h " + type(term).__name__)
        ST["steps"] += 1
        n = ST["steps"]
        if n == 1:
            sync()
            res["setup_s"] = time.perf_counter() - T0
        if n == A.skip + 1:
            sync()
            c.update(t=time.perf_counter(), dm=ST["dm"], tries=ST["tries"])
            c["t0"], c["dm0"], c["tries0"] = c["t"], c["dm"], c["tries"]
            ST["phase"] = "plain"
        elif n == A.skip + A.plain // 2 + 1:
            sync()
            c["half"] = {"wall_s": time.perf_counter() - c["t0"], "stages": ST["dm"] - c["dm0"],
                         "tries": ST["tries"] - c["tries0"], "sample_intervals": A.plain // 2}
            c.update(t1=time.perf_counter(), dm1=ST["dm"], tries1=ST["tries"])
        elif n == A.skip + A.plain + 1:
            sync()
            res["plain"] = {"wall_s": time.perf_counter() - c["t"], "stages": ST["dm"] - c["dm"],
                            "tries": ST["tries"] - c["tries"], "sample_intervals": A.plain}
            res["plain"]["halves"] = [c["half"], {"wall_s": time.perf_counter() - c["t1"], "stages": ST["dm"] - c["dm1"],
                                                  "tries": ST["tries"] - c["tries1"],
                                                  "sample_intervals": A.plain - A.plain // 2}]
            if torch.cuda.is_available():
                res["peak_alloc_GiB"] = torch.cuda.max_memory_allocated() / 2**30
                res["peak_reserved_GiB"] = torch.cuda.max_memory_reserved() / 2**30
            else:
                res["peak_alloc_GiB"] = res["peak_reserved_GiB"] = float("nan")
            if A.window == 0:
                report_plain(res)
                raise SystemExit(0)
            c["prof"] = profile(activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA])
            c["prof"].__enter__()
            c.update(t=time.perf_counter(), dm=ST["dm"], tries=ST["tries"])
            ST["phase"] = "prof"
        if ST["phase"] == "prof":
            with record_function("L:LLGSolver.step"):
                r = f_step(self, state, dt, **k)
        else:
            r = f_step(self, state, dt, **k)
        if A.trace:
            Ms = state.material["Ms"]
            res.setdefault("trace", []).append((float(state.t), ST["tries"], float((Ms * state.m[..., 0:1]).mean()),
                                                int(state._step)))
        if A.window and n == A.skip + A.plain + A.window:          # end of the profiled window (none with --window 0)
            sync()
            wall = time.perf_counter() - c["t"]
            n_dm, n_try = ST["dm"] - c["dm"], ST["tries"] - c["tries"]
            c["prof"].__exit__(None, None, None)
            ST["phase"] = "done"
            res["window"] = analyse(c["prof"], wall, n_dm, n_try, A.window)
            report(res)
            raise SystemExit(0)
        return r

    magnumnp.LLGSolver.step = step
    run_loops.main(rl)
    print("[prof_step] the run ended before the end of the window (%d sample intervals)" % ST["steps"], flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
