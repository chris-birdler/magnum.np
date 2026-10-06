#!/usr/bin/env python
"""
Post-processing on the instance during the Sobol design (PLAN 4.0; standing go of Chris 2026-10-06).
Loop: for every finished run (DONE) with snapshots (*.vti) and without POSTPROC_DONE:
  1. phasemap_batch.py and phase_extra.py (phasemap.json, phase_extra.json, slices.npz)
  2. check: both JSON files exist and hold the same number of stages as there are snapshot stages
  3. md5 of the three result files -> postproc.md5 (in the runs folder), then delete the *.vti of this run
     (ONLY *.vti: init.pt and checkpoint.pt are never touched here)
  4. touch POSTPROC_DONE
A run whose check fails keeps its snapshots and is listed in postproc_errors.log.
Stops when the file --stop exists and no run is left to process.

    python instance_postproc.py --runs runs --stop runs/ALL10_DONE --workers 4
"""
import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def process(d):
    stages = {re.match(r"m_(.+)_c\d+_ph\d+\.vti", f.name).group(1) for f in d.glob("m_*_c*_ph*.vti")}
    for script in ("phasemap_batch.py", "phase_extra.py"):
        r = subprocess.run([sys.executable, str(HERE / script), str(d)], capture_output=True, text=True,
                           env=dict(os.environ, OMP_NUM_THREADS="6"))
        if r.returncode != 0:
            return d, "%s failed: %s" % (script, r.stderr[-500:])
    try:
        n1 = len(json.loads((d / "phasemap.json").read_text())["stages"])
        n2 = len(json.loads((d / "phase_extra.json").read_text())["stages"])
    except Exception as e:                                  # noqa: BLE001
        return d, "result files not readable: %s" % e
    if not (n1 == n2 == len(stages)) or not (d / "slices.npz").exists():
        return d, "stage count mismatch: snapshots %d, phasemap %d, phase_extra %d" % (len(stages), n1, n2)
    lines = ["%s  %s/%s\n" % (md5(d / f), d.name, f) for f in ("phasemap.json", "phase_extra.json", "slices.npz")]
    return d, lines


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--stop", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--once", action="store_true", help="one pass, then exit (test)")
    a = ap.parse_args()
    runs = pathlib.Path(a.runs)
    with cf.ProcessPoolExecutor(a.workers) as ex:
        while True:
            todo = [d for d in sorted(runs.iterdir()) if d.is_dir() and (d / "DONE").exists()
                    and not (d / "POSTPROC_DONE").exists() and any(d.glob("m_*.vti"))]
            for d, res in ex.map(process, todo):
                if isinstance(res, str):
                    with open(runs / "postproc_errors.log", "a") as f:
                        f.write("%s %s: %s\n" % (time.strftime("%F %T"), d.name, res))
                    (d / "POSTPROC_FAILED").write_text(res)
                    continue
                with open(runs / "postproc.md5", "a") as f:
                    f.writelines(res)
                for v in d.glob("*.vti"):
                    v.unlink()
                (d / "POSTPROC_DONE").write_text(time.strftime("%F %T\n"))
                print("%s %s: analysed, snapshots deleted" % (time.strftime("%T"), d.name), flush=True)
            if a.once or (pathlib.Path(a.stop).exists() and not todo):
                break
            time.sleep(60)


if __name__ == "__main__":
    main()
