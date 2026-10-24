#!/usr/bin/env python
"""
Verdict of D4 (PLAN 14b): outlier test 3 on the patched tree. Reads the json files of glitch_count --reps --json
(<stem>_redo_<file>.json) of the P arm and of the rate control, and the logs of the P arm (also of the fillers that were
stopped by their PID).

- P arm: each complete run must have 0 tries with a difference and 0 non-finite values; the logs must have no
  "redo event" line.
- Rate control: n_C = events above the thresholds of outlier test 3 whose first part (glitch_count --term_rel 1e-9) is
  "spmv"; exposure E_C = try evaluations of the complete runs (H2_d53 pooled).
- Test: p = (E_C / (E_C + E_P))^n_C, pass at p <= 0.01 (with the P arm free of differences).

    python d4_eval.py --p P1.json P2.json ... --c C1.json C2.json ... --logs H3_d53P*.log [--json out.json]
"""
import argparse
import json
import pathlib
import re

THR3 = (1e-6, 1e-12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", nargs="+", required=True, help="glitch_count redo json files of the complete P runs")
    ap.add_argument("--c", nargs="+", required=True, help="glitch_count redo json files of the complete rate-control runs")
    ap.add_argument("--logs", nargs="*", default=[], help="all logs of the P arm, stopped fillers included")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--json")
    a = ap.parse_args()
    P = [json.loads(pathlib.Path(f).read_text()) for f in a.p]
    C = [json.loads(pathlib.Path(f).read_text()) for f in a.c]
    E_P = sum(d["try_evaluations"] for d in P)
    p_diff = sum(d["tries_with_difference"] for d in P)
    p_bad = sum(d["non_finite"] for d in P)
    log_ev = {f: len(re.findall(r"redo event", pathlib.Path(f).read_text())) for f in a.logs}
    E_C = sum(d["try_evaluations"] for d in C)
    n_C, other = 0, []
    for d in C:
        for ev in d["event_list"]:
            above = not (ev["derr"] <= THR3[0] and ev["dx1"] <= THR3[1])
            first = (ev.get("first") or {}).get("term")
            if above and first == "spmv":
                n_C += 1
            elif above:
                other.append((d["file"], ev["n"], first))
    p = (E_C / (E_C + E_P)) ** n_C if E_C + E_P > 0 else float("nan")
    clean = p_diff == 0 and p_bad == 0 and sum(log_ev.values()) == 0
    verdict = "PASS" if clean and p <= a.alpha else ("FAIL (difference in the P arm)" if not clean
                                                     else "NOT DETERMINED (p above %g)" % a.alpha)
    for d in P:
        print("P    %-12s try evaluations %7d, tries with a difference %d, non-finite %d" % (
            d["file"], d["try_evaluations"], d["tries_with_difference"], d["non_finite"]))
    for f, k in log_ev.items():
        print("P    log %-40s redo event lines %d" % (pathlib.Path(f).name, k))
    for d in C:
        n = sum(1 for ev in d["event_list"] if not (ev["derr"] <= THR3[0] and ev["dx1"] <= THR3[1])
                and (ev.get("first") or {}).get("term") == "spmv")
        print("rate %-12s try evaluations %7d, events above the thresholds of test 3 with spmv first %d" % (
            d["file"], d["try_evaluations"], n))
    for o in other:
        print("rate event above the thresholds with another first part: %s n %d first %s" % o)
    print("E_P %d, E_C %d, n_C %d, p = (E_C / (E_C + E_P))^n_C = %.3e -> %s" % (E_P, E_C, n_C, p, verdict))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps({"E_P": E_P, "E_C": E_C, "n_C": n_C, "p": p, "verdict": verdict,
                                                    "p_tries_with_difference": p_diff, "p_non_finite": p_bad,
                                                    "p_log_events": log_ev, "rate_other_events": other}, indent=1))


if __name__ == "__main__":
    main()
