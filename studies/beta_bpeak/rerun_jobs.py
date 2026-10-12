#!/usr/bin/env python
"""
Job lines of PLAN 4.5 from the original job lines.

  --reruns   every finished run of the given job files with a stage whose calm tail is < 2 (calm cycle:
             w_loop > 0 and |w_dis / w_loop - 1| <= 0.20; measured cycles as in analyze.py) gets a rerun X_mc2:
             the original arguments (with the original --init_from, if any; without one, the rerun relaxes in the
             same process as the original did: the same start path gives the same cycles, PROTOKOLL 7.59), --min_calm 2,
             --max_cycles_per_amp replaced (FILE:MAX), all --snap_* removed. With --cap_vh: jobs in file order
             until the estimated cost (V100-h per cycle of the original × MAX cycles × stages) reaches the cap;
             the others are listed, not written.
  --cond     conditioning test (PLAN 4.5 D): for each given run X three runs cond_X_b9 / _b20 / _b35 with the
             first 1 / 2 / 3 amplitudes of X, --cycles_last_amp 30, --min_calm 2, --max_cycles_per_amp 15,
             no --init_from (C0 … C5 relaxed in the process; PROTOKOLL 7.59), no snapshots
  --retire   before a queue: a folder X_mc2 that started from init.pt while its job line in the queue file has no
             --init_from is renamed to X_mc2init (kept as a sensitivity record, not used by the analysis)

    python rerun_jobs.py --reruns jobs/mesh300.txt:12 jobs/meshK3.txt:12 jobs/sobol_main.txt:15 \
        jobs/sobol_corners.txt:15 --out jobs/rerun_mc2.txt
    python rerun_jobs.py --cond C0 C1 C2 C3 C4 C5 --source jobs/sobol_main.txt --out jobs/cond.txt
"""
import argparse
import json
import pathlib
import shlex

from analyze import steady_tail                   # light import (no torch); = run_loops.calm_tail, test

HERE = pathlib.Path(__file__).resolve().parent
MIN_CALM = 2
CALM_TOL = 0.20


def job_lines(path):
    out = []
    for line in open(path):
        if line.strip() and not line.startswith("#"):
            name, _, args = line.strip().partition(" ")
            out.append((name, args))
    return out


def set_arg(tok, key, values):
    """Replace (or add) --key with values (list of strings; [] = flag); None removes it."""
    out, i = [], 0
    while i < len(tok):
        if tok[i] == "--" + key:
            i += 1
            while i < len(tok) and not tok[i].startswith("--"):
                i += 1
            continue
        out.append(tok[i])
        i += 1
    if values is not None:
        out += ["--" + key] + list(values)
    return out


def get_arg(tok, key):
    if "--" + key not in tok:
        return None
    i = tok.index("--" + key) + 1
    vals = []
    while i < len(tok) and not tok[i].startswith("--"):
        vals.append(tok[i])
        i += 1
    return vals


def strip_snap(tok):
    for k in ("snap_phases", "snap_cycles", "snap_mT"):
        tok = set_arg(tok, k, None)
    return tok


def short_stages(d):
    """[(stage index, calm tail)] of the measured stages of a finished run with a calm tail < MIN_CALM."""
    st = json.loads((d / "summary.json").read_text())["stages"]
    out = []
    for k, s in enumerate(x for x in st if x.get("measure")):
        meas = s["cycles"][2:] if "h_amp_initial" in s else s["cycles"][1:]
        ct = len(steady_tail(meas, CALM_TOL))
        if ct < MIN_CALM:
            out.append((k, ct))
    return out


def h_per_cycle(d):
    st = json.loads((d / "summary.json").read_text())["stages"]
    c = [x for s in st for x in s["cycles"]]
    return sum(x["wall_s"] for x in c) / 3600.0 / max(len(c), 1)


def reruns(specs, runs, cap_vh):
    lines, skipped, cost = [], [], 0.0
    for spec in specs:
        path, mx = spec.rsplit(":", 1)
        for name, args in job_lines(path):
            d = runs / name
            if not (d / "DONE").exists():
                continue
            short = short_stages(d)
            if not short:
                continue
            tok = shlex.split(args)
            n_st = len(get_arg(tok, "b_list") or ["x"])
            est = h_per_cycle(d) * int(mx) * n_st
            if cap_vh is not None and cost + est > cap_vh:
                skipped.append("%s (estimate %.1f V100-h; stages with calm tail < 2: %s)" % (name, est, short))
                continue
            cost += est
            tok = strip_snap(tok)
            tok = set_arg(tok, "max_cycles_per_amp", [mx])
            tok = set_arg(tok, "min_calm", [str(MIN_CALM)])
            lines.append("%s_mc2 %s" % (name, " ".join(tok)))
    return lines, skipped, cost


def cond(names, source, amps=(9, 20, 35), max_prev=15):
    """Conditioning runs: the design history up to each tested amplitude, then 30 cycles at it. max_prev caps the
    earlier stages (15 at α 0.1, PLAN 4.5; 9 = the design length at α 0.01, decision Chris 2026-10-12)."""
    src = dict(job_lines(source))
    lines = []
    for n in names:
        tok = shlex.split(src[n])
        bl = get_arg(tok, "b_list")
        for k, mT in enumerate(amps):
            t = strip_snap(list(tok))
            t = set_arg(t, "b_list", bl[:k + 1])
            t = set_arg(t, "n_amp", [str(k + 1)])
            t = set_arg(t, "cycles_last_amp", ["30"])
            t = set_arg(t, "max_cycles_per_amp", [str(max_prev)])
            t = set_arg(t, "min_calm", [str(MIN_CALM)])
            lines.append("cond_%s_b%d %s" % (n, mT, " ".join(t)))
    return lines


def retire(queue_file, runs):
    """Rename X_mc2 -> X_mc2init (…init2, …) when the folder started from init.pt but its queue line has none."""
    out = []
    for name, args in job_lines(queue_file):
        d = runs / name
        if not name.endswith("_mc2") or not (d / "config.json").exists():
            continue
        started_from_init = json.loads((d / "config.json").read_text()).get("init_from") is not None
        if started_from_init and get_arg(shlex.split(args), "init_from") is None:
            k, new = 1, runs / (name + "init")
            while new.exists():
                k += 1
                new = runs / ("%sinit%d" % (name, k))
            d.rename(new)
            out.append("%s -> %s" % (name, new.name))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--reruns", nargs="*", default=None, help="FILE:MAX_CYCLES ...")
    ap.add_argument("--cap_vh", type=float, default=None)
    ap.add_argument("--cond", nargs="*", default=None, help="run names for the conditioning test")
    ap.add_argument("--source", default=None, help="job file of the --cond runs")
    ap.add_argument("--cond_amps", nargs="*", type=int, default=[9, 20, 35], help="tested amplitudes (mT)")
    ap.add_argument("--cond_max_prev", type=int, default=15, help="cycle cap of the earlier stages")
    ap.add_argument("--retire", default=None, help="queue file: rename the init.pt variants (see above), then stop")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.retire:
        moved = retire(a.retire, pathlib.Path(a.runs))
        print("retired: %s" % (", ".join(moved) or "none"))
        return
    lines = []
    if a.reruns:
        lines, skipped, cost = reruns(a.reruns, pathlib.Path(a.runs), a.cap_vh)
        print("%d reruns, upper bound %.1f V100-h (all stages at the maximum)%s" % (len(lines), cost, ("; not written (cap): " + "; ".join(skipped))
                                                      if skipped else ""))
    if a.cond:
        lines += cond(a.cond, a.source, tuple(a.cond_amps), a.cond_max_prev)
        print("%d conditioning runs" % (len(a.cond_amps) * len(a.cond)))
    pathlib.Path(a.out).write_text("".join(l + "\n" for l in lines))
    print("written %s (%d lines)" % (a.out, len(lines)))


if __name__ == "__main__":
    main()
