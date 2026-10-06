#!/usr/bin/env python
"""
Post-processing on the instance during the Sobol design (PLAN 4.0; standing go of Chris 2026-10-06).
Loop: for every finished run (DONE) with snapshots (*.vti) or accumulators (acc_*.npz, run_loops --acc) and
without POSTPROC_DONE (accumulators: acc_reduce.py -> acc.json, acc_slices.npz; checked: one stage per file, sum
check p_box / p_dis = 1, R1 finite; then the acc_*.npz are deleted like the snapshots, decision Chris 2026-10-06):
  1. phasemap_batch.py, phase_extra.py and domains_batch.py (phasemap.json, phase_extra.json, slices.npz,
     domains.json: wall measures, switched volume, l_C of the last cycle per stage, of init.pt and of the
     end state checkpoint.pt)
  2. check: the three JSON files exist and hold the same number of stages as there are snapshot stages
  3. md5 of the four result files -> postproc.md5 (in the runs folder), then delete the *.vti of this run
     (ONLY *.vti: init.pt and checkpoint.pt are never touched here)
  4. touch POSTPROC_DONE
A run whose check fails keeps its snapshots and is listed in postproc_errors.log.
Only runs whose name matches --pattern (default: the Sobol runs S###, S###_a001, C#, C#_a001, K#_#) are
touched. Every error of one run is caught (logged, the run keeps its snapshots); the loop never stops on it.
Stops when the file --stop exists and no run is left to process.

    python instance_postproc.py --runs runs --stop runs/ALL10_DONE --workers 4
"""
import argparse
import concurrent.futures as cf
from concurrent.futures.process import BrokenProcessPool
import hashlib
import json
import math
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
    """Analyse one finished run: snapshots (*.vti) and / or accumulators (acc_*.npz). Returns (d, md5 lines) or
    (d, error text)."""
    stages = {re.match(r"m_(.+)_c\d+_ph\d+\.vti", f.name).group(1) for f in d.glob("m_*_c*_ph*.vti")}
    accs = [f for f in d.glob("acc_*_*.npz") if not f.name.endswith(".tmp.npz")]
    scripts = (["phasemap_batch.py", "phase_extra.py", "domains_batch.py"] if stages else []) + \
        (["acc_reduce.py"] if accs else [])
    for script in scripts:
        r = subprocess.run([sys.executable, str(HERE / script), str(d)], capture_output=True, text=True,
                           env=dict(os.environ, OMP_NUM_THREADS="6"))
        if r.returncode != 0:
            return d, "%s failed: %s" % (script, r.stderr[-500:])
    files = []
    if stages:
        try:
            pm = json.loads((d / "phasemap.json").read_text())["stages"]
            n1 = len(pm)
            n2 = len(json.loads((d / "phase_extra.json").read_text())["stages"])
            n3 = len(json.loads((d / "domains.json").read_text())["stages"])
        except Exception as e:                                  # noqa: BLE001
            return d, "result files not readable: %s" % e
        bad = [s.get("stage") for s in pm if not all(math.isfinite(s.get(k, float("nan"))) for k in ("R1", "R_np", "p_dis"))]
        if bad:
            return d, "phasemap.json: R1 / R_np / p_dis not finite in stage(s) %s" % bad
        if not (n1 == n2 == n3 == len(stages)) or not (d / "slices.npz").exists():
            return d, "stage count mismatch: snapshots %d, phasemap %d, phase_extra %d, domains %d" % (
                len(stages), n1, n2, n3)
        files += ["phasemap.json", "phase_extra.json", "slices.npz", "domains.json"]
    if accs:
        try:
            ac = json.loads((d / "acc.json").read_text())["stages"]
        except Exception as e:                                  # noqa: BLE001
            return d, "acc.json not readable: %s" % e
        if len(ac) != len(accs) or not (d / "acc_slices.npz").exists():
            return d, "accumulators: %d files, acc.json %d stages" % (len(accs), len(ac))
        bad = [s["stage"] for s in ac if not (s["check_ok"] and math.isfinite(s["R1"]))]
        if bad:
            return d, "accumulators: sum check failed or R1 not finite in stage(s) %s" % bad
        files += ["acc.json", "acc_slices.npz"]
    return d, ["%s  %s/%s\n" % (md5(d / f), d.name, f) for f in files]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--stop", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--once", action="store_true", help="one pass, then exit (test)")
    ap.add_argument("--pattern", default=r"^(S\d{3}(_a001)?|C\d(_a001)?|K\d_\d)$",
                    help="regular expression of the run names that may be processed (and whose *.vti are deleted)")
    a = ap.parse_args()
    runs = pathlib.Path(a.runs)
    pat = re.compile(a.pattern)
    with cf.ProcessPoolExecutor(a.workers) as ex:
        while True:
            todo = [d for d in sorted(runs.iterdir()) if d.is_dir() and pat.match(d.name) and (d / "DONE").exists()
                    and not (d / "POSTPROC_DONE").exists() and not (d / "POSTPROC_FAILED").exists()
                    and (any(d.glob("m_*.vti")) or any(d.glob("acc_*_*.npz")))]
            futs = {ex.submit(process, d): d for d in todo}
            for fu in cf.as_completed(futs):
                d = futs[fu]
                try:
                    d, res = fu.result()
                except BrokenProcessPool:
                    raise                                   # a worker was killed (e.g. OOM): not the run's fault;
                    #                                         the supervisor loop of the chain restarts this script
                except Exception as e:                      # noqa: BLE001 - one run must not stop the loop
                    res = "exception: %r" % e
                if isinstance(res, str):
                    with open(runs / "postproc_errors.log", "a") as f:
                        f.write("%s %s: %s\n" % (time.strftime("%F %T"), d.name, res))
                    (d / "POSTPROC_FAILED").write_text(res)
                    continue
                try:
                    with open(runs / "postproc.md5", "a") as f:
                        f.writelines(res)
                    for v in list(d.glob("*.vti")) + list(d.glob("acc_*_*.npz")):
                        v.unlink()
                    (d / "POSTPROC_DONE").write_text(time.strftime("%F %T\n"))
                    print("%s %s: analysed, snapshots deleted" % (time.strftime("%T"), d.name), flush=True)
                except Exception as e:                      # noqa: BLE001
                    with open(runs / "postproc_errors.log", "a") as f:
                        f.write("%s %s: after analysis: %r\n" % (time.strftime("%F %T"), d.name, e))
            if a.once or (pathlib.Path(a.stop).exists() and not todo):
                break
            time.sleep(60)


if __name__ == "__main__":
    main()
