"""Shared helpers for the magnum.np benchmark scripts."""
import json
import os
import resource
import subprocess
import sys
import time

import torch


def set_dtype(dtype_arg):
    """Select f32/f64, falling back to the raw torch default for magnum.np
    versions that do not provide set_precision (so baselines can be benchmarked)."""
    import magnumnp as mag
    if dtype_arg == "f32":
        if hasattr(mag, "set_precision"):
            mag.set_precision("single")
        else:
            torch.set_default_dtype(torch.float32)


def time_fn(fn, warmup=5, iters=50):
    """Return average wall time per call in ms (cuda-event based on GPU)."""
    for _ in range(warmup):
        fn()
    if torch.cuda.is_available():
        torch.cuda.synchronize()
        start = torch.cuda.Event(enable_timing=True)
        end = torch.cuda.Event(enable_timing=True)
        start.record()
        for _ in range(iters):
            fn()
        end.record()
        torch.cuda.synchronize()
        return start.elapsed_time(end) / iters
    t0 = time.perf_counter()
    for _ in range(iters):
        fn()
    return (time.perf_counter() - t0) / iters * 1000.


def reset_peak_memory():
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()


def peak_memory():
    mem = {"ru_maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024., 1)}
    if torch.cuda.is_available():
        mem["cuda_peak_mb"] = round(torch.cuda.max_memory_allocated() / 1e6, 1)
    return mem


def run_variants(script, variants, out, name):
    """Run `script --variant <args>` in a fresh subprocess per variant and
    collect the json line each one prints; write the records to <out>/<name>.json."""
    records = []
    for variant in variants:
        cmd = [sys.executable, script, "--variant"] + [str(v) for v in variant]
        print(">>", " ".join(cmd), flush=True)
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr)
            raise RuntimeError("variant %s failed" % (variant,))
        record = json.loads(result.stdout.strip().splitlines()[-1])
        print("  ", json.dumps(record), flush=True)
        records.append(record)

    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, name + ".json"), "w") as fd:
        json.dump(records, fd, indent=2)
    return records


def out_dir_from_args():
    args = sys.argv[1:]
    if len(args) >= 2 and args[0] == "--out":
        return args[1]
    print("usage: %s --out <results-dir>" % sys.argv[0])
    sys.exit(1)
