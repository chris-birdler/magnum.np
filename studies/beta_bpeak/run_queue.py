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
"""
import argparse
import os
import pathlib
import queue
import shlex
import subprocess
import sys
import threading
import time

HERE = pathlib.Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobfile")
    ap.add_argument("--gpus", default="0", help="comma separated CUDA device ids")
    ap.add_argument("--runs", default=str(HERE / "runs"), help="output root")
    ap.add_argument("--shard", default="0/1", help="i/n: run every n-th job starting at i")
    ap.add_argument("--dry", action="store_true", help="print the commands only")
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
            out.mkdir(parents=True, exist_ok=True)
            env = dict(os.environ, CUDA_DEVICE=str(gpu))
            t0 = time.time()
            print("[start] gpu %s  %s" % (gpu, name), flush=True)
            with open(out / "stdout.log", "a") as log:
                rc = subprocess.call(cmd, env=env, stdout=log, stderr=subprocess.STDOUT)
            print("[%s] gpu %s  %s  %.2f h" % ("done" if rc == 0 else "FAIL rc=%d" % rc, gpu, name,
                                                (time.time() - t0) / 3600), flush=True)

    threads = [threading.Thread(target=worker, args=(g.strip(),)) for g in a.gpus.split(",")]
    for t in threads:
        t.start()
    for t in threads:
        t.join()


if __name__ == "__main__":
    main()
