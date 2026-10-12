"""Task 1: inventory and integrity of every local run folder (audit 2026-10-12).

Output: audit_20261012/inventory.csv (one line per run folder), inventory.md (tables per group).
Checks per run: DONE, relaxation converged (summary.json init), stages complete, job line found and equal to the
config, stdout.log equal to summary.json (last printed value of every cycle), calm tail per stage (rule 7.54,
>= 2), rerun replacement and determinism of the rerun against its original (cycle by cycle up to the first
extended stage), git commit.
"""
import collections
import csv
import json
import math
import pathlib
import re

import numpy as np

import audit_lib as al

RUNS = al.RUNS
JOBDIRS = [al.STUDY / "jobs", RUNS / "mc2b_logs"]


def group_of(n):
    rules = [
        (r"^S\d{3}$", "design: Sobol main", "7.44/7.54/7.65"),
        (r"^C\d$", "design: centre replicate", "7.44/7.65"),
        (r"^K\d_\d$", "design: corner replicate", "7.44/7.65"),
        (r"^(S\d{3}|C\d|K\d_\d)_mc2$", "design: rerun to 2 calm cycles", "7.65"),
        (r"^C\d_a001$", "damping test alpha 0.01 (centre)", "7.48/7.62"),
        (r"^cond_C\d_b\d+$", "conditioning test (PLAN 4.5 D)", "7.65"),
        (r"^M300_dx(3|15)_s\d+_b\d+$", "mesh300 (d/l_ex 300, base)", "7.25/7.64/7.65"),
        (r"^M300_dx(3|15)_s\d+_b\d+_mc2$", "mesh300 rerun", "7.65"),
        (r"^MK3_dx(3|15)_s\d+_b(9|50|100)$", "K3 mesh check (9/50/100 mT)", "7.52/7.65"),
        (r"^MK3_dx(3|15)_s\d+_b(20|35|70)$", "K3 fill-in (20/35/70 mT)", "7.57/7.65"),
        (r"^MK3_dx(3|15)_s\d+_b\d+_mc2$", "K3 rerun", "7.65"),
        (r"^MK3cg_dx3_s\d+_b\d+$", "K3 same-state test", "7.63"),
        (r"^M212_", "mesh check d/l_ex 212 (PLAN 4.6)", "pending"),
        (r"^B300_s\d+$", "baseline d/l_ex 300 (base point)", "7.17"),
        (r"^PH_s\d+$", "phase map / ring-down", "7.26/7.27"),
        (r"^F9_(dbl|half)_s\d+$", "frequency test 9/20 mT (f/2, 2f)", "7.26"),
        (r"^V_s\d+$", "wall snapshots", "7.18"),
        (r"^A010c_s\d+$", "alpha 0.01 clean (base)", "7.17/7.22"),
        (r"^X(010|100)_s\d+$", "relaxation test", "7.15"),
        (r"^A0(10|30)_s\d+$", "alpha test (old relaxation)", "7.12"),
        (r"^F_(base|corner)_(f|2f)_s\d+$", "frequency test (old relaxation)", "7.9"),
        (r"^R300f(01)?_s\d+_b\d+$", "surface: fraction model", "7.29"),
        (r"^R300v01_s\d+_b\d+$", "paired alpha 0.01 (base)", "7.31"),
        (r"^R212_dx(3|15)_s\d+$", "mesh d/l_ex 212 (old relaxation)", "7.10"),
        (r"^S(212|300)_dx3f_s\d+(_b\d+)?$", "surface roughness (old relaxation)", "7.11"),
        (r"^P1_", "planning pilot block 1", "7.8"),
        (r"^P2_", "planning pilot block 2", "7.8"),
        (r"^P3_|^R_dx15_", "planning pilot block 3 (mesh d/l_ex 96)", "7.8"),
        (r"^M_(L12|floor)_dx", "mesh test d/l_ex 96", "7.1"),
        (r"^S_(L12|floor|L12_a010)_dx", "steady test d/l_ex 96", "7.2"),
        (r"^F_(do|si)_atol", "numerical floor test", "7.3"),
        (r"^I_", "initial-state test", "7.4"),
        (r"^D_", "demagnetisation test", "7.5"),
        (r"^C_dx", "mesh convergence d/l_ex 96", "7.6"),
        (r"^bench", "benchmark", "-"),
    ]
    for pat, g, sec in rules:
        if re.match(pat, n):
            return g, sec
    return "unknown", "-"


def parse_jobs():
    jobs = collections.defaultdict(list)
    for jd in JOBDIRS:
        for f in sorted(jd.glob("*.txt")):
            for line in f.read_text().splitlines():
                t = line.split()
                if not t or t[0].startswith("#"):
                    continue
                args, k = {}, None
                for tok in t[1:]:
                    if tok.startswith("--"):
                        k = tok[2:]
                        args[k] = []
                    elif k is not None:
                        args[k].append(tok)
                jobs[t[0]].append((f.name, args))
    return jobs


def cmp_args(args, cfg):
    bad = []
    for k, v in args.items():
        if k in ("job_timeout_h", "init_wait_h", "no_resume", "allow_new_commit", "gpu"):
            continue
        if k not in cfg:
            bad.append("%s:not in config" % k)
            continue
        c = cfg[k]
        if not v:                                    # flag
            if c is not True:
                bad.append("%s:%r" % (k, c))
            continue
        if k == "init_from":
            if c is None or not str(c).endswith(v[0]):
                bad.append("init_from:%r vs %s" % (c, v[0]))
            continue
        vals = v if isinstance(c, list) else v[:1]
        cl = c if isinstance(c, list) else [c]
        if len(vals) != len(cl):
            bad.append("%s:len" % k)
            continue
        for a, b in zip(vals, cl):
            try:
                fa = float(a)
                if b is None or abs(fa - float(b)) > 1e-9 * max(1.0, abs(fa)):
                    bad.append("%s:%s vs %r" % (k, a, b))
            except ValueError:
                if str(a) != str(b):
                    bad.append("%s:%s vs %r" % (k, a, b))
    return bad


LINE = re.compile(r"^\[(\w+) c(\d+)\].*?w_loop=([-\d.eE+]+)\s+w_dis=([-\d.eE+]+)")


def stdout_check(run, summ):
    f = run / "stdout.log"
    if not f.exists():
        return "no stdout"
    last = {}
    for line in f.read_text(errors="replace").splitlines():
        m = LINE.match(line.strip())
        if m:
            last[(m.group(1), int(m.group(2)))] = (float(m.group(3)), float(m.group(4)))
    n, bad = 0, 0
    for st in summ["stages"]:
        for c in st["cycles"]:
            k = (st["name"], c["cycle"])
            if k not in last:
                continue
            n += 1
            wl, wd = last[k]
            if abs(wl - c["w_loop"]) > 6e-4 * abs(c["w_loop"]) + 1e-30 or abs(wd - c["w_dis"]) > 6e-4 * abs(c["w_dis"]) + 1e-30:
                bad += 1
    tot = sum(len(st["cycles"]) for st in summ["stages"])
    return "ok %d/%d" % (n, tot) if bad == 0 else "MISMATCH %d of %d" % (bad, n)


def gate(orig, rerun):
    """cycle-by-cycle equality of rerun and original up to the first stage with another cycle count."""
    _, _, so = al.read(orig)
    _, _, sr = al.read(rerun)
    n, maxrel = 0, 0.0
    for a, b in zip(so["stages"], sr["stages"]):
        for ca, cb in zip(a["cycles"], b["cycles"]):
            for key in ("w_loop", "w_dis", "b_peak"):
                maxrel = max(maxrel, abs(ca[key] - cb[key]) / max(abs(ca[key]), 1e-300))
            n += 1
        if len(a["cycles"]) != len(b["cycles"]):
            break
    return n, maxrel


def main():
    jobs = parse_jobs()
    rows = []
    for d in sorted(RUNS.iterdir()):
        if not d.is_dir() or d.name in ("mc2b_logs",):
            continue
        if d.name.startswith("_archive"):
            n_arch = sum(1 for x in d.iterdir() if x.is_dir())
            rows.append({"name": d.name, "group": "archive: Sobol v2 old stress axes", "section": "7.36/7.38",
                         "status": "archived (old non-isotropic stress axes, not used)", "notes": "%d run folders" % n_arch})
            continue
        g, sec = group_of(d.name)
        rec = {"name": d.name, "group": g, "section": sec}
        try:
            run, cfg, summ = al.read(d)
        except Exception as e:                        # noqa: BLE001
            rec.update(status="unreadable", notes=str(e))
            rows.append(rec)
            continue
        rec["done"] = (d / "DONE").exists()
        init = (summ or {}).get("init", {})
        rec["relax_mode"] = cfg.get("relax_mode")
        rec["relax_conv"] = init.get("relax_converged")
        rec["relax_it"] = init.get("relax_iterations")
        rec["E_relax"] = init.get("E_over_Kd")
        rec["n60_relax"] = init.get("n_pairs_gt60")
        rec["init_from"] = cfg.get("init_from")
        for k in ("d_lex", "dx_lex", "Leff_lex", "r_p", "phi", "alpha", "seed"):
            rec[k] = cfg.get(k)
        rec["f_MHz"] = round(cfg.get("freq", float("nan")) / 1e6, 3)
        rec["commit"] = (cfg.get("git_commit") or "")[:8]
        rec["gpu"] = cfg.get("gpu")
        if summ is None:
            rec.update(status="no summary", notes="")
            rows.append(rec)
            continue
        st_meas = [st for st in summ["stages"] if st.get("measure")]
        nb = len(cfg.get("b_list") or []) or cfg.get("n_amp")
        rec["n_stages"] = len(st_meas)
        rec["stages_complete"] = all(st.get("complete", False) and st["cycles"] for st in st_meas) and len(st_meas) >= (nb or 0)
        tails, evs = [], 0
        for st in st_meas:
            cyc = al.measured(st)
            t = al.tail(cyc)
            tails.append(len(t))
            evs += sum(al.is_event(c) for c in cyc)
        rec["min_calm_tail"] = min(tails) if tails else None
        rec["n_stage_tail_lt2"] = sum(t < 2 for t in tails)
        rec["n_event_cycles"] = evs
        rec["n_meas_cycles"] = sum(len(al.measured(st)) for st in st_meas)
        rec["stdout"] = stdout_check(d, summ)
        jl = jobs.get(d.name, [])
        rec["jobfile"] = ";".join(sorted({f for f, _ in jl}))
        if jl:
            mism = [cmp_args(a, cfg) for _, a in jl]
            best = min(mism, key=len)
            rec["job_vs_config"] = "ok" if not best else "DIFF " + ",".join(best)
        else:
            rec["job_vs_config"] = "no job line"
        # status
        notes = []
        if d.name.endswith("_mc2"):
            orig = d.name[:-4]
            n, mr = gate(orig, d.name)
            notes.append("gate vs %s: %d cycles, max rel diff %.1e" % (orig, n, mr))
        if g.startswith("design") or g.startswith("damping") or g.startswith("conditioning") or \
                g.startswith("mesh300") or g.startswith("K3"):
            if not rec["done"]:
                status = "INCOMPLETE"
            elif not rec["relax_conv"]:
                status = "INVALID: relaxation not converged"
            elif al.rerun_of(d.name):
                status = "replaced by rerun %s" % al.rerun_of(d.name)
            else:
                status = "current, valid"
                if rec["n_stage_tail_lt2"] and not g.startswith("damping"):
                    status = "current, %d stage(s) with calm tail < 2" % rec["n_stage_tail_lt2"]
        elif g in ("baseline d/l_ex 300 (base point)", "phase map / ring-down", "frequency test 9/20 mT (f/2, 2f)",
                   "wall snapshots", "alpha 0.01 clean (base)", "relaxation test", "surface: fraction model",
                   "paired alpha 0.01 (base)"):
            if not rec["done"]:
                status = "stopped / incomplete (not used)"
            elif rec["relax_conv"]:
                status = "valid (converged), earlier study; steady rule not applied in its section"
            else:
                status = "INVALID: relaxation not converged"
        elif g == "benchmark":
            status = "test (benchmark, not used)"
        else:
            status = ("withdrawn (7.13: relaxation not converged)" if not rec["relax_conv"]
                      else "test / pilot")
        rec["status"] = status
        rec["notes"] = "; ".join(notes)
        rows.append(rec)
    keys = ["name", "group", "section", "status", "done", "relax_mode", "relax_conv", "relax_it", "E_relax", "n60_relax",
            "init_from", "d_lex", "dx_lex", "Leff_lex", "r_p", "phi", "alpha", "seed", "f_MHz", "n_stages",
            "stages_complete", "n_meas_cycles", "n_event_cycles", "min_calm_tail", "n_stage_tail_lt2", "stdout",
            "jobfile", "job_vs_config", "commit", "gpu", "notes"]
    with open(al.HERE / "inventory.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})
    # summary per group
    by = collections.defaultdict(list)
    for r in rows:
        by[(r["group"], r.get("section", ""))].append(r)
    L = ["# Inventory of the local run folders (audit 2026-10-12)", "",
         "%d entries in runs/ (folders only; files and logs not counted)." % len(rows), "",
         "| group | PROTOKOLL | n | status counts | relax converged | DONE | stdout = summary | job line = config | min calm tail |",
         "|---|---|---|---|---|---|---|---|---|"]
    for (g, sec), rs in sorted(by.items(), key=lambda kv: kv[0]):
        sc = collections.Counter(r.get("status", "") for r in rs)
        conv = sum(1 for r in rs if r.get("relax_conv"))
        done = sum(1 for r in rs if r.get("done"))
        so = sum(1 for r in rs if str(r.get("stdout", "")).startswith("ok"))
        jc = sum(1 for r in rs if r.get("job_vs_config") == "ok")
        mct = [r.get("min_calm_tail") for r in rs if isinstance(r.get("min_calm_tail"), int)]
        L.append("| %s | %s | %d | %s | %d | %d | %d | %d | %s |" % (
            g, sec, len(rs), "; ".join("%s: %d" % kv for kv in sc.items()), conv, done, so, jc,
            min(mct) if mct else "-"))
    L += ["", "## Problems", ""]
    for r in rows:
        probs = []
        if str(r.get("stdout", "")).startswith("MISMATCH"):
            probs.append("stdout " + r["stdout"])
        if str(r.get("job_vs_config", "")).startswith("DIFF"):
            probs.append(r["job_vs_config"])
        if r.get("job_vs_config") == "no job line" and not r["group"].startswith(("archive", "benchmark")):
            probs.append("no job line")
        if r.get("stages_complete") is False:
            probs.append("stages incomplete")
        if probs:
            L.append("- %s (%s, %s): %s" % (r["name"], r["group"], r.get("status"), "; ".join(probs)))
    L += ["", "## Reruns (gate recomputed locally)", ""]
    for r in rows:
        if r["name"].endswith("_mc2"):
            L.append("- %s: %s; init_from %s" % (r["name"], r["notes"], r.get("init_from")))
    (al.HERE / "inventory.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
