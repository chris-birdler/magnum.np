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
  minus the bytes already in the run folder. Thus the running jobs cannot fill the disk; the queue waits instead.
"""
import argparse
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
    return cells * (n_snap * B_SNAP_CELL + 2 * B_PT_CELL)


def folder_bytes(p):
    return sum(f.stat().st_size for f in p.glob("*") if f.is_file()) if p.exists() else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobfile")
    ap.add_argument("--gpus", default="0", help="comma separated CUDA device ids")
    ap.add_argument("--runs", default=str(HERE / "runs"), help="output root")
    ap.add_argument("--shard", default="0/1", help="i/n: run every n-th job starting at i")
    ap.add_argument("--dry", action="store_true", help="print the commands only")
    ap.add_argument("--min_free_gb", type=float, default=0.0,
                    help="disk guard with reservation (see above): free space left after all reservations (0 = off)")
    a = ap.parse_args()

    i_sh, n_sh = (int(x) for x in a.shard.split("/"))
    jobs = []
    for k, line in enumerate(l for l in open(a.jobfile) if l.strip() and not l.startswith("#")):
        if k % n_sh != i_sh:
            continue
        name, _, args = line.strip().partition(" ")
        jobs.append((name, args))

    runs = pathlib.Path(a.runs)
    q = queue.Queue()
    for name, args in jobs:
        if (runs / name / "DONE").exists():
            print("[skip] %s (DONE)" % name)
            continue
        q.put((name, args))

    running = {}                                    # name -> (run folder, expected bytes)
    lock = threading.Lock()

    def reserve(name, args, out):
        """Disk guard with reservation; returns when the job may start (and registers it)."""
        need = expected_bytes(args)
        waited = False
        while True:
            with lock:
                pending = sum(max(e - folder_bytes(o), 0) for o, e in running.values())
                free = shutil.disk_usage(str(runs)).free
                if free - pending - need >= a.min_free_gb * 1e9:
                    running[name] = (out, need)
                    if waited:
                        print("[go] %s: free %.1f GB, reserved %.1f GB + %.1f GB" % (name, free / 1e9, pending / 1e9,
                                                                                   need / 1e9), flush=True)
                    return
            if not waited:
                print("[wait] %s: free %.1f GB - reserved %.1f GB - needs %.1f GB < %.0f GB" %
                      (name, free / 1e9, pending / 1e9, need / 1e9, a.min_free_gb), flush=True)
                waited = True
            time.sleep(300)

    def worker(gpu):
        while True:
            try:
                name, args = q.get_nowait()
            except queue.Empty:
                return
            out = runs / name
            cmd = [sys.executable, str(HERE / "run_loops.py")] + shlex.split(args) + ["--out", str(out)]
            if a.dry:
                print("CUDA_DEVICE=%s %s" % (gpu, " ".join(cmd)))
                continue
            if a.min_free_gb > 0:
                reserve(name, args, out)
            out.mkdir(parents=True, exist_ok=True)
            env = dict(os.environ, CUDA_DEVICE=str(gpu))
            t0 = time.time()
            print("[start] gpu %s  %s" % (gpu, name), flush=True)
            with open(out / "stdout.log", "a") as log:
                rc = subprocess.call(cmd, env=env, stdout=log, stderr=subprocess.STDOUT)
            with lock:
                running.pop(name, None)
            print("[%s] gpu %s  %s  %.2f h" % ("done" if rc == 0 else "FAIL rc=%d" % rc, gpu, name,
                                                (time.time() - t0) / 3600), flush=True)

    threads = [threading.Thread(target=worker, args=(g.strip(),)) for g in a.gpus.split(",")]
    for t in threads:
        t.start()
    for t in threads:
        t.join()


if __name__ == "__main__":
    main()
