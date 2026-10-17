import json, math, pathlib, re, sys
import numpy as np
sys.path.insert(0, "/tmp/claude-1000/-home-christoph-Documents-code-magnum-np-powder/f96c028a-8b43-4c77-ac5c-7ee7423ab0ef/scratchpad")
from scan_rules import RUNS, stage_cycles, make_rule, eps, RULES  # noqa

base = make_rule()


def main_dirs():
    pat = re.compile(r"^(S\d{3}|C\d|K\d_\d)$")
    out = []
    for d in sorted(RUNS.iterdir()):
        if d.is_dir() and pat.match(d.name):
            r = d.parent / (d.name + "_mc2")
            d = r if (r / "DONE").exists() else d
            if json.loads((d / "config.json").read_text())["alpha"] == 0.1:
                out.append(d)
    return out


def tail_measures(dirs, label):
    """Per stage over the baseline tail: within-tail dW, tail-mean eps, SD of ln w, drift of ln w (slope x (n-1))."""
    rows = {k: [] for k in range(7)}
    for d in dirs:
        s, sts = stage_cycles(d)
        for k, (st, cyc) in enumerate(sts):
            t = base(cyc, None)
            if not t:
                continue
            w = np.array([c["w_loop"] for c in t]); wd = np.array([c["w_dis"] for c in t])
            lw = np.log(w)
            dw_in = np.abs(np.diff(w)) / w[1:] if len(w) > 1 else np.array([])
            em = wd.sum() / w.sum() - 1.0
            sd = lw.std(ddof=1) if len(w) > 1 else np.nan
            drift = np.polyfit(np.arange(len(w)), lw, 1)[0] * (len(w) - 1) if len(w) > 2 else np.nan
            rows[min(k, 6)].append(dict(n=len(w), dwin=dw_in, em=em, sd=sd, drift=drift))
    print("\n## %s: per-stage tail measures (baseline tail)" % label)
    print("amp nst  dW_in(50/90/99/max)       |eps_mean| 50/90/max  n>0.02 n>0.05   SD lnw 50/90/max  n>0.05   |drift| 50/90/max n>0.05")
    tot = dict(e2=0, e5=0, s5=0, d5=0, n=0)
    for k in range(7):
        r = rows[k]
        if not r:
            continue
        dwin = np.concatenate([x["dwin"] for x in r]) if any(len(x["dwin"]) for x in r) else np.array([np.nan])
        em = np.abs(np.array([x["em"] for x in r])); sd = np.array([x["sd"] for x in r]); dr = np.abs(np.array([x["drift"] for x in r]))
        q = lambda x, p=(0.5, 0.9, 0.99, 1.0): "/".join("%.2g" % v for v in np.nanquantile(x, p))
        print("%3d %3d  %-24s  %-18s %4d %4d   %-17s %4d   %-17s %4d" % (
            k, len(r), q(dwin), q(em, (0.5, 0.9, 1.0)), (em > 0.02).sum(), (em > 0.05).sum(),
            q(sd, (0.5, 0.9, 1.0)), np.nansum(sd > 0.05), q(dr, (0.5, 0.9, 1.0)), np.nansum(dr > 0.05)))
        tot["e2"] += (em > 0.02).sum(); tot["e5"] += (em > 0.05).sum(); tot["s5"] += np.nansum(sd > 0.05); tot["d5"] += np.nansum(dr > 0.05); tot["n"] += len(r)
    print("total stages %d: |eps_mean|>0.02 %d, >0.05 %d; SD lnw>0.05 %d; |drift|>0.05 %d" % (tot["n"], tot["e2"], tot["e5"], tot["s5"], tot["d5"]))


def dlnw_by_amp(dirs, label):
    """Per stage: ln(mean w of alt tail) - ln(mean w of baseline tail), by amplitude; only stages with >=1 alt cycle."""
    print("\n## %s: per-stage d ln w (alt tail - baseline tail), mean +- SE (n stages with a change / n stages), by amplitude" % label)
    for name, kw in RULES.items():
        if name == "eps0.20":
            continue
        alt = make_rule(**kw)
        per = {k: [] for k in range(7)}
        for d in dirs:
            s, sts = stage_cycles(d)
            for k, (st, cyc) in enumerate(sts):
                t0 = base(cyc, None); t1 = alt(cyc, None)
                if t0 and t1:
                    per[min(k, 6)].append(math.log(np.mean([c["w_loop"] for c in t1]) / np.mean([c["w_loop"] for c in t0])))
        parts = []
        for k in range(7):
            x = np.array(per[k])
            parts.append("%+.4f+-%.4f(%d/%d)" % (x.mean(), x.std(ddof=1) / math.sqrt(len(x)), int((np.abs(x) > 1e-12).sum()), len(x)))
        print("%-8s " % name + " ".join(parts))


def new_runs():
    print("\n## NEW-SETTING runs: per measured cycle (first_measured on), eps / closure / dW_rel / b/b_target; steps, wall")
    groups = [("TB", sorted(RUNS.glob("TB_*"))), ("N1g", sorted(RUNS.glob("N1g_*"))),
              ("N0b tight", sorted(RUNS.glob("N0b_F_s94[12]_atol*"))), ("N0b stage", sorted(d for d in RUNS.glob("N0b_*") if re.match(r"N0b_(F|I3|I15|T3|T15)_s94\d$", d.name)))]
    for g, ds in groups:
        print("-- %s" % g)
        for d in ds:
            cfg = json.loads((d / "config.json").read_text())
            if not (d / "summary.json").exists():
                print("%-28s no summary" % d.name); continue
            s, sts = stage_cycles(d)
            tag = "dx %.2g atol %s/%s rtol %s sf %s" % (cfg["dx_lex"], cfg.get("atol"), cfg.get("atol_cycles"), cfg.get("rtol_cycles"), cfg.get("strayfield"))
            for st, cyc in sts:
                allc = st["cycles"]
                s1 = " ".join("c%d:e%+.4f/cl%.1e/dW%.3f" % (c["cycle"], eps(c), c["closure"], c.get("dW_rel", np.nan)) for c in cyc)
                bt = np.mean([c["b_peak"] for c in cyc]) / st["b_target"] if cyc else np.nan
                print("%-30s %s | ncyc %d calm_tail %s | b/bt %.3f | steps/cyc %d wall/cyc %.0f s | %s" % (
                    d.name, tag, len(allc), st.get("calm_tail"), bt, np.mean([c["steps"] for c in allc]), np.mean([c["wall_s"] for c in allc]), s1))
            ini = s.get("init", {})
            h = ini.get("relax_history", [])
            if h:
                print("    relax: conv %s it %s last rel_dE %s  torque rms/max %.3g/%.3g  E/Kd %.6g" % (
                    ini.get("relax_converged"), ini.get("relax_iterations"), ["%.1e" % x["rel_dE"] for x in h[-3:]], h[-1]["torque_rms"], h[-1]["torque_max"], h[-1]["E_over_Kd"]))


if __name__ == "__main__":
    m = main_dirs()
    tail_measures(m, "OLD MAIN")
    l6 = sorted(d for d in RUNS.iterdir() if re.match(r"^K1_\d_L6$", d.name))
    tail_measures(l6, "K1_k_L6 (atol 1e-6)")
    a001 = sorted(d for d in RUNS.iterdir() if re.match(r"^C\d_a001$", d.name))
    tail_measures(a001, "ALPHA 0.01")
    dlnw_by_amp(m, "OLD MAIN")
    new_runs()
