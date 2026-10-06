#!/usr/bin/env python
"""
Run a job file on several GPUs of one machine (one job per GPU at a time).

    python run_queue.py jobs/production.txt --gpus 0,1,2,3
    python run_queue.py jobs/production.txt --gpus 0,1,2,3,4,5,6,7 --shard 0/2   # half of the jobs
    python run_queue.py jobs/bench.txt --gpus 0 --dry

* A job is skipped if runs/<name>/DONE exists.
* An interrupted job restarts from its checkpoint when you start the queue
  again (run_loops.py resumes automatically).
* stdout/stderr of each job: runs/<name>/stdout.log
* --shard i/n takes every n-th job starting at i (to split one job file over
  several vast.ai instances without a shared file system).
* --min_free_gb G (disk guard with reservation, PLAN 4.0): a job starts only if
      free disk - (expected bytes still to come of all running jobs) - (expected bytes of the new job) >= G.
  Expected bytes of a job (`expected_bytes`): snapshots (phases x cycles x amplitudes x N^3 x 15 B, measured
  47.7 MB for N = 148 in single precision) + init.pt + checkpoint.pt (2 x N^3 x 12 B). Still to come = expected
  minus the bytes already in the run folder. With --acc: + amplitudes x N^3 x 40 B (accumulators). Thus the running jobs cannot fill the disk; the queue waits instead.
Unattended operation (audit 3, item 4):
* a failed job (exit code != 0, timeout, exception) goes once more to the end of the queue (--retries 1); it
  resumes from its checkpoint. Exit code 3 (relaxation gate not met) is deterministic (same seed, same
  iterations): no plain retry, but once again with --relax_maxiter --gate_maxiter (200000) and --no_resume
  (audit 3, item 5); if it fails again the job is reported as excluded.
* failure burst: failures of 2 DIFFERENT jobs on one GPU within 10 min (a job that fails twice is the job's
  fault, not the GPU's) pause this GPU for 30 min; after the 3rd burst the GPU takes no more jobs (defective
  GPU). Each burst goes to <runs>/ALERT; the other GPUs go on.
* --job_timeout_h H: a job that runs longer is killed (FAIL "timeout").
* every exception of a worker is caught; the reservation of the job is always released.
* --status FILE: a JSON status (running, done, failed, waiting, disabled GPUs, disk, reservations, alerts,
  post-processing state) is written every 5 min and at the end.
* exit code 0 only if all jobs are DONE.
* the retry state of each job (attempts, gate retry with its args, excluded) is kept in
  <runs>/<name>/QUEUE_STATE.json, so a restart of the queue (instance restart) does not give a job new
  attempts (audit 3, C2): an excluded job is skipped; a job in its gate retry resumes with the raised
  relax_maxiter (without --no_resume). Delete the file by hand to give a job new attempts.
"""
import argparse
import json
import os
import pathlib
import queue
import shlex
import shutil
import subprocess
import sys
import threading
import time

HERE = pathlib.Path(__file__).resolve().parent
B_SNAP_CELL = 15.0      # bytes per cell of one snapshot (vti, fp32, measured 47.7 MB / 148^3 on the instance)
B_PT_CELL = 12.0        # bytes per cell of init.pt / checkpoint.pt (measured 38.9 MB / 148^3)
B_ACC_CELL = 40.0       # bytes per cell of one acc_*.npz (p 4 + m0 12 + m1 24; particle cells only: upper bound)


def job_args(args):
    """Option dict of a job line: --key v1 v2 ... -> {key: [v1, v2, ...]}."""
    d, key = {}, None
    for t in shlex.split(args):
        if t.startswith("--"):
            key = t[2:]
            d[key] = []
        elif key is not None:
            d[key].append(t)
    return d


def expected_bytes(args):
    """Upper estimate of the bytes that a job writes into its run folder (snapshots + state files)."""
    from anisotropy_model import box_edge
    d = job_args(args)
    a_lex, _ = box_edge(float(d["d_lex"][0]), float(d["phi"][0]))
    cells = round(a_lex / float(d["dx_lex"][0])) ** 3
    phases = int(d.get("snap_phases", ["0"])[0])
    n_snap = 0
    if phases > 0:
        n_amp = len(d["snap_mT"]) if "snap_mT" in d else len(d.get("b_list", [])) or int(d.get("n_amp", ["7"])[0])
        n_snap = phases * int(d.get("snap_cycles", ["1"])[0]) * n_amp
    n_acc = (len(d.get("b_list", [])) or int(d.get("n_amp", ["7"])[0])) if "acc" in d else 0
    return cells * (n_snap * B_SNAP_CELL + 2 * B_PT_CELL + n_acc * B_ACC_CELL)


def folder_bytes(p):
    """Bytes of the files in a run folder. A file can vanish while it is counted (atomic save: x.tmp -> x;
    post-processing deletes snapshots): such a file is skipped."""
    n = 0
    for f in (p.glob("*") if p.exists() else []):
        try:
            n += f.stat().st_size
        except OSError:
            pass
    return n


def default_command(args, out):
    return [sys.executable, str(HERE / "run_loops.py")] + shlex.split(args) + ["--out", str(out)]


def postproc_state(runs):
    st = {"done": 0, "failed": [], "waiting": 0}
    for d in runs.iterdir() if runs.exists() else []:
        if not d.is_dir():
            continue
        if (d / "POSTPROC_DONE").exists():
            st["done"] += 1
        elif (d / "POSTPROC_FAILED").exists():
            st["failed"].append(d.name)
        elif (d / "DONE").exists() and any(d.glob("m_*.vti")):
            st["waiting"] += 1
    return st


def load_queue_state(d):
    try:
        return json.loads((d / "QUEUE_STATE.json").read_text())
    except (OSError, ValueError):
        return {}


def save_queue_state(d, st):
    d.mkdir(parents=True, exist_ok=True)
    tmp = d / "QUEUE_STATE.json.tmp"
    tmp.write_text(json.dumps(st, indent=1))
    os.replace(tmp, d / "QUEUE_STATE.json")


def gate_args(args, maxiter):
    """Job args for the retry after a gate failure: --relax_maxiter raised, --no_resume (the unconverged
    init.pt of the first attempt must not be reused)."""
    t = shlex.split(args)
    if "--relax_maxiter" in t:
        t[t.index("--relax_maxiter") + 1] = str(maxiter)
    else:
        t += ["--relax_maxiter", str(maxiter)]
    if "--no_resume" not in t:
        t.append("--no_resume")
    return " ".join(shlex.quote(x) for x in t)


def run_jobs(jobs, gpus, runs, min_free_gb=0.0, retries=1, job_timeout_h=0.0, status=None, command=default_command,
             burst_n=2, burst_s=600.0, pause_s=1800.0, burst_max=3, status_s=300.0, poll_s=300.0, gate_maxiter=200000):
    """Run the (name, args) jobs on the GPUs; returns {"done": [...], "failed": {name: reason}}."""
    runs = pathlib.Path(runs)
    runs.mkdir(parents=True, exist_ok=True)
    q = queue.Queue()
    attempts, failed, done, disabled, alerts, gate_retried = {}, {}, [], set(), [], set()
    for name, args in jobs:
        if (runs / name / "DONE").exists():
            print("[skip] %s (DONE)" % name, flush=True)
            continue
        qs = load_queue_state(runs / name)
        if qs.get("excluded"):
            failed[name] = "excluded earlier: %s" % qs["excluded"]
            print("[skip] %s (%s)" % (name, failed[name]), flush=True)
            continue
        attempts[name] = int(qs.get("attempts", 0))
        if qs.get("gate_retried"):
            gate_retried.add(name)
            args = " ".join(shlex.quote(x) for x in shlex.split(qs["args"]) if x != "--no_resume")
        q.put((name, args))
    lock = threading.Lock()
    running = {}                                    # name -> (run folder, expected bytes, gpu, start time)
    active = {"n": 0}                               # jobs taken from the queue and not finished
    t_start = time.time()

    def alert(msg):
        line = "%s %s" % (time.strftime("%F %T"), msg)
        alerts.append(line)
        print("[ALERT] " + msg, flush=True)
        with open(runs / "ALERT", "a") as f:
            f.write(line + "\n")

    def write_status(final=False):
        if not status:
            return
        with lock:
            free = shutil.disk_usage(str(runs)).free
            pending = sum(max(e - folder_bytes(o), 0) for o, e, _, _ in running.values())
            st = {"time": time.strftime("%F %T"), "final": final, "hours": (time.time() - t_start) / 3600,
                  "running": {n: {"gpu": g, "hours": (time.time() - t0) / 3600} for n, (_, _, g, t0) in running.items()},
                  "done": len(done), "failed": dict(failed), "waiting": q.qsize(), "disabled_gpus": sorted(disabled),
                  "free_gb": free / 1e9, "reserved_gb": pending / 1e9, "alerts": list(alerts),
                  "postproc": postproc_state(runs)}
        tmp = pathlib.Path(str(status) + ".tmp")
        tmp.write_text(json.dumps(st, indent=1))
        os.replace(tmp, status)

    def reserve(name, args, out, gpu):
        need = expected_bytes(args) if min_free_gb > 0 else 0
        waited = False
        while True:
            with lock:
                pending = sum(max(e - folder_bytes(o), 0) for o, e, _, _ in running.values())
                free = shutil.disk_usage(str(runs)).free
                if min_free_gb <= 0 or free - pending - need >= min_free_gb * 1e9:
                    running[name] = (out, need, gpu, time.time())
                    if waited:
                        print("[go] %s: free %.1f GB, reserved %.1f GB + %.1f GB" % (name, free / 1e9, pending / 1e9,
                                                                                   need / 1e9), flush=True)
                    return
            if not waited:
                print("[wait] %s: free %.1f GB - reserved %.1f GB - needs %.1f GB < %.0f GB" %
                      (name, free / 1e9, pending / 1e9, need / 1e9, min_free_gb), flush=True)
                waited = True
            time.sleep(poll_s)

    def run_one(name, args, gpu):
        out = runs / name
        reserve(name, args, out, gpu)
        try:
            out.mkdir(parents=True, exist_ok=True)
            env = dict(os.environ, CUDA_DEVICE=str(gpu))
            print("[start] gpu %s  %s" % (gpu, name), flush=True)
            with open(out / "stdout.log", "a") as log:
                proc = subprocess.Popen(command(args, out), env=env, stdout=log, stderr=subprocess.STDOUT)
                try:
                    rc = proc.wait(timeout=job_timeout_h * 3600 if job_timeout_h > 0 else None)
                    return rc, "exit code %d" % rc
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                    return -1, "timeout %.1f h" % job_timeout_h
        except Exception as e:                      # noqa: BLE001 - a worker must not die
            return -2, "exception %r" % e
        finally:
            with lock:
                running.pop(name, None)

    def worker(gpu):
        fails, bursts = [], 0                       # fails: (time, job name) on this GPU
        while True:
            try:
                name, args = q.get(timeout=5)
            except queue.Empty:
                with lock:
                    if active["n"] == 0:            # nothing running that could come back to the queue
                        return
                continue
            with lock:
                active["n"] += 1
            t0 = time.time()
            rc, why = run_one(name, args, gpu)
            h = (time.time() - t0) / 3600
            with lock:
                if rc == 0 and (runs / name / "DONE").exists():
                    done.append(name)
                    failed.pop(name, None)
                    print("[done] gpu %s  %s  %.2f h" % (gpu, name, h), flush=True)
                else:
                    failed[name] = why
                    if rc == 3:
                        retry = gate_maxiter > 0 and name not in gate_retried
                        if retry:
                            gate_retried.add(name)
                            args = gate_args(args, gate_maxiter)
                            failed[name] = "relaxation gate; again with relax_maxiter %d" % gate_maxiter
                        else:
                            failed[name] = "relaxation gate%s: excluded" % (" (also with relax_maxiter %d)"
                                                                            % gate_maxiter if gate_maxiter > 0 else "")
                    else:
                        attempts[name] = attempts.get(name, 0) + 1
                        retry = attempts[name] <= retries
                    print("[FAIL] gpu %s  %s  %.2f h  %s%s" % (gpu, name, h, why, ", again later" if retry else ""),
                          flush=True)
                    save_queue_state(runs / name, {"attempts": attempts.get(name, 0),
                                                    "gate_retried": name in gate_retried, "args": args,
                                                    "excluded": None if retry else failed[name]})
                    if retry:
                        q.put((name, args))
                    fails = [(t, n) for t, n in fails if time.time() - t < burst_s] + [(time.time(), name)]
                active["n"] -= 1
            if len({n for _, n in fails}) >= burst_n:
                bursts += 1
                fails = []
                if bursts >= burst_max:
                    with lock:
                        disabled.add(gpu)
                    alert("gpu %s disabled: burst %d of failures of %d different jobs within %.0f min (last: %s %s)"
                          % (gpu, bursts, burst_n, burst_s / 60, name, why))
                    if len(disabled) == len(gpus):
                        alert("all GPUs disabled: the queue stops")
                    return
                alert("gpu %s paused %.0f min: burst %d of failures of %d different jobs within %.0f min (last: %s %s)"
                      % (gpu, pause_s / 60, bursts, burst_n, burst_s / 60, name, why))
                time.sleep(pause_s)

    stop = threading.Event()

    def status_loop():
        while not stop.wait(status_s):
            try:
                write_status()
            except Exception as e:                  # noqa: BLE001
                print("[status] %r" % e, flush=True)

    threads = [threading.Thread(target=worker, args=(g,)) for g in gpus]
    st_thread = threading.Thread(target=status_loop, daemon=True)
    st_thread.start()
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    stop.set()
    left = []
    while not q.empty():
        left.append(q.get()[0])
    for n in left:
        failed.setdefault(n, "not run (no GPU left)")
    write_status(final=True)
    return {"done": done, "failed": {n: w for n, w in failed.items() if n not in done}}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobfile")
    ap.add_argument("--gpus", default="0", help="comma separated CUDA device ids")
    ap.add_argument("--runs", default=str(HERE / "runs"), help="output root")
    ap.add_argument("--shard", default="0/1", help="i/n: run every n-th job starting at i")
    ap.add_argument("--dry", action="store_true", help="print the commands only")
    ap.add_argument("--min_free_gb", type=float, default=0.0,
                    help="disk guard with reservation (see above): free space left after all reservations (0 = off)")
    ap.add_argument("--retries", type=int, default=1, help="extra attempts of a failed job (not for exit code 3)")
    ap.add_argument("--job_timeout_h", type=float, default=0.0, help="kill a job after this time (0 = off)")
    ap.add_argument("--status", default=None, help="JSON status file, written every 5 min")
    ap.add_argument("--gate_maxiter", type=int, default=200000,
                    help="relax_maxiter of the one retry after a gate failure (exit code 3); 0 = no retry")
    a = ap.parse_args()

    i_sh, n_sh = (int(x) for x in a.shard.split("/"))
    jobs = []
    for k, line in enumerate(l for l in open(a.jobfile) if l.strip() and not l.startswith("#")):
        if k % n_sh != i_sh:
            continue
        name, _, args = line.strip().partition(" ")
        jobs.append((name, args))
    gpus = [g.strip() for g in a.gpus.split(",")]
    if a.dry:
        for name, args in jobs:
            print("CUDA_DEVICE=%s %s" % (gpus[0], " ".join(default_command(args, pathlib.Path(a.runs) / name))))
        return
    res = run_jobs(jobs, gpus, a.runs, a.min_free_gb, a.retries, a.job_timeout_h, a.status,
                   gate_maxiter=a.gate_maxiter)
    print("[end] done %d, failed %d: %s" % (len(res["done"]), len(res["failed"]),
                                           ", ".join("%s (%s)" % kv for kv in sorted(res["failed"].items()))), flush=True)
    sys.exit(0 if not res["failed"] else 1)


if __name__ == "__main__":
    main()
