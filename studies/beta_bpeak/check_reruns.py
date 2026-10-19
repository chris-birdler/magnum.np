#!/usr/bin/env python
"""
Determinism gate of PLAN 4.5: a rerun must repeat its reference run cycle by cycle until the protocols part.
  X_mc2          (rerun with --min_calm 2)       reference X
  cond_X_bNN     (conditioning test, PLAN 4.5 D) reference X
Stage by stage, the first min(n_ref, n_new) cycles are compared (w_loop, w_dis, b_peak; relative difference
<= RTOL). The first stage where the numbers of cycles differ is the last stage compared: after it the start states
differ by design. A stage of the rerun that is still running is compared as far as it goes.

    python check_reruns.py [--runs runs]                    # one check, exit code 1 on a mismatch
    python check_reruns.py --watch 600 --stop_file runs/STOP_NEW_JOBS --until runs/MC2_QUEUE_DONE
                                                            # every 600 s; a mismatch creates the stop file
"""
import argparse
import json
import pathlib
import re
import time

RTOL = 1e-6
KEYS = ("w_loop", "w_dis", "b_peak")


def reference(name):
    m = re.match(r"^(.+)_mc2$", name) or re.match(r"^cond_(.+)_b\d+$", name)
    return m.group(1) if m else None


def stages(d):
    """The finished stages of a run (summary.json is written at the end of each stage)."""
    try:
        return json.loads((d / "summary.json").read_text())["stages"]
    except (OSError, ValueError):
        return []


def compare(ref, new):
    """-> (n_cycles_compared, mismatch text or None)."""
    n = 0
    for k, (so, sn) in enumerate(zip(ref, new)):
        co, cn = so["cycles"], sn["cycles"]
        for i in range(min(len(co), len(cn))):
            for key in KEYS:
                a, b = co[i][key], cn[i][key]
                if abs(a - b) > RTOL * max(abs(a), abs(b), 1e-300):
                    return n, "stage %d (%s) cycle %d: %s %.9g vs %.9g" % (k, so["name"], i, key, a, b)
            n += 1
        if len(co) != len(cn):
            break
    return n, None


def check(runs):
    out, bad = [], []
    for d in sorted(runs.iterdir()):
        ref = reference(d.name) if d.is_dir() else None
        if ref is None or not (runs / ref / "summary.json").exists():
            continue
        n, why = compare(stages(runs / ref), stages(d))
        out.append("%s vs %s: %d cycles equal%s" % (d.name, ref, n, "" if why is None else "; MISMATCH " + why))
        if why is not None:
            bad.append(d.name)
    return out, bad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=str(pathlib.Path(__file__).resolve().parent / "runs"))
    ap.add_argument("--watch", type=float, default=0.0, help="repeat every this many seconds")
    ap.add_argument("--stop_file", default=None)
    ap.add_argument("--until", default=None, help="end the watch after a last check when this file exists")
    a = ap.parse_args()
    runs = pathlib.Path(a.runs)
    while True:
        last = a.until is not None and pathlib.Path(a.until).exists()
        lines, bad = check(runs)
        stamp = time.strftime("%F %T")
        if bad:
            print("%s check_reruns: MISMATCH in %s" % (stamp, ", ".join(bad)), flush=True)
            for line in lines:
                print("  " + line, flush=True)
            if a.stop_file:
                pathlib.Path(a.stop_file).write_text("check_reruns %s: %s\n" % (stamp, ", ".join(bad)))
            raise SystemExit(1)
        if not a.watch or last:
            print("%s check_reruns: %d runs, all equal" % (stamp, len(lines)), flush=True)
            for line in lines:
                print("  " + line, flush=True)
            return
        time.sleep(a.watch)


if __name__ == "__main__":
    main()
