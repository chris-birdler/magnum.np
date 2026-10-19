"""Review scan: per-cycle eps, closure, dW_rel in the kept cycles; alternative steady rules; effect on ln w and beta.
Read-only on the run folders."""
import json, math, pathlib, sys
import numpy as np

CODE = pathlib.Path("/home/christoph/Documents/code/magnum.np-newell/studies/beta_bpeak")
RUNS = pathlib.Path("/home/christoph/Documents/code/magnum.np-beta_bpeak/studies/beta_bpeak/runs")
sys.path.insert(0, str(CODE))
import analyze, analyze_design as ad  # noqa

AMPS = ad.AMPS_MT


def eps(c):
    return c["w_dis"] / c["w_loop"] - 1.0 if c["w_loop"] != 0 else float("inf")


def make_rule(eps_tol=0.20, clo_tol=None, dw_tol=None):
    def viol(c):
        if c["w_loop"] <= 0 or abs(eps(c)) > eps_tol:
            return True
        if clo_tol is not None and c["closure"] > clo_tol:
            return True
        if dw_tol is not None and c.get("dW_rel", 0.0) > dw_tol:
            return True
        return False

    def tail(cycles, tol):
        ev = [i for i, c in enumerate(cycles) if viol(c)]
        return cycles[ev[-1] + 1:] if ev else cycles
    return tail


RULES = {"eps0.20": dict(), "eps0.10": dict(eps_tol=0.10), "eps0.05": dict(eps_tol=0.05), "eps0.02": dict(eps_tol=0.02),
         "clo1e-3": dict(clo_tol=1e-3), "clo3e-4": dict(clo_tol=3e-4), "clo1e-4": dict(clo_tol=1e-4),
         "dW0.10": dict(dw_tol=0.10), "dW0.05": dict(dw_tol=0.05)}


def stage_cycles(d):
    s = json.loads((d / "summary.json").read_text())
    out = []
    for st in s["stages"]:
        if not st.get("measure") or not st.get("cycles"):
            continue
        fm = analyze.first_measured(st) if "h_amp_initial" in st else analyze.SKIP
        out.append((st, st["cycles"][fm:]))
    return s, out


def kept_stats(dirs, label):
    """eps, closure, dW_rel of the cycles kept by the baseline rule, per nominal amplitude index."""
    base = make_rule()
    per = {k: {"eps": [], "clo": [], "dW": []} for k in range(7)}
    for d in dirs:
        s, sts = stage_cycles(d)
        for k, (st, cyc) in enumerate(sts):
            for c in base(cyc, None):
                per[k]["eps"].append(eps(c)); per[k]["clo"].append(c["closure"]); per[k]["dW"].append(c.get("dW_rel", np.nan))
    print("\n## %s: kept cycles (baseline |eps|<=0.2 tail), quantiles 50/90/99/max" % label)
    print("amp  n   eps(50/90/99/max, signed median)          closure(50/90/99/max)            dW_rel(50/90/99/max)")
    for k in range(7):
        e = np.array(per[k]["eps"]); cl = np.array(per[k]["clo"]); dw = np.array(per[k]["dW"])
        if len(e) == 0:
            continue
        q = lambda x: "/".join("%.2g" % v for v in np.nanquantile(np.abs(x), [0.5, 0.9, 0.99, 1.0]))
        print("%3d %4d  |eps| %s med %+.4f   clo %s   dW %s" % (k, len(e), q(e), np.median(e), q(cl), q(dw)))
    return per


def rule_counts(dirs, label):
    print("\n## %s: alternative rules (tail after last violating cycle)" % label)
    print("rule      dropped/total cycles   stages with 0 kept / <2 kept   runs with a missing stage (0 kept)")
    for name, kw in RULES.items():
        tail = make_rule(**kw)
        tot = drop = z = lt2 = 0
        runs0 = set()
        for d in dirs:
            s, sts = stage_cycles(d)
            for st, cyc in sts:
                t = tail(cyc, None)
                tot += len(cyc); drop += len(cyc) - len(t)
                z += len(t) == 0; lt2 += len(t) < 2
                if len(t) == 0:
                    runs0.add(d.name)
        print("%-9s %5d / %5d (%4.1f %%)       %4d / %4d                    %d" % (name, drop, tot, 100 * drop / max(tot, 1), z, lt2, len(runs0)))


def run_outputs(dirs, rule_kw, steady_min=1):
    analyze.STEADY_TAIL = 0.2
    analyze.STEADY_MIN = steady_min
    analyze.steady_tail = make_rule(**rule_kw)
    recs = {}
    for d in dirs:
        try:
            r = ad.run_record(d)
        except Exception as e:  # noqa
            r = None
        if r is not None:
            recs[d.name] = r
    return recs


def ols_range(recs, key):
    names = [n for n in recs if np.isfinite(recs[n][key])]
    X = ad.coded([recs[n] for n in names])
    y = np.array([recs[n][key] for n in names])
    T = ad.terms(X)
    b, *_ = np.linalg.lstsq(T, y, rcond=None)
    eff = []
    for j in range(4):
        Xp, Xm = X.copy(), X.copy(); Xp[:, j] = 1; Xm[:, j] = -1
        eff.append(float(np.mean(ad.terms(Xp) @ b - ad.terms(Xm) @ b)))
    return eff, len(names)


def compare_rules(dirs, label, keys=("beta_all", "lnw@10mT", "lnw@50mT", "lnw@100mT", "beta@9-35mT", "beta@70-150mT"), design=True):
    base = run_outputs(dirs, {})
    print("\n## %s: per-run change against the baseline rule (mean +- SE over runs, max |d|, n runs with a change)" % label)
    for name, kw in RULES.items():
        if name == "eps0.20":
            continue
        for smin in (1, 2):
            alt = run_outputs(dirs, kw, smin)
            parts = []
            for k in keys:
                d = np.array([alt[n][k] - base[n][k] for n in base if n in alt and np.isfinite(alt[n][k]) and np.isfinite(base[n][k])])
                nan_new = sum(1 for n in base if n in alt and np.isfinite(base[n][k]) and not np.isfinite(alt[n][k]))
                nz = int(np.sum(np.abs(d) > 1e-12))
                parts.append("%s %+.4f+-%.4f max %.3f nchg %d nlost %d" % (k, d.mean(), d.std(ddof=1) / math.sqrt(len(d)), np.abs(d).max(), nz, nan_new))
            print("%-8s min%d | " % (name, smin) + " | ".join(parts))
            if design and smin == 1:
                for k in ("beta_all", "lnw@50mT"):
                    e0, n0 = ols_range(base, k); e1, n1 = ols_range(alt, k)
                    print("          OLS15 range effect %s (lnQ, r_p, phi, lnd): base %s (n %d)  alt %s (n %d)" % (
                        k, " ".join("%+.3f" % v for v in e0), n0, " ".join("%+.3f" % v for v in e1), n1))
    analyze.steady_tail = make_rule()


if __name__ == "__main__":
    import re
    pat = re.compile(r"^(S\d{3}|C\d|K\d_\d)$")
    orig = sorted([d for d in RUNS.iterdir() if d.is_dir() and pat.match(d.name)])
    main = []
    for d in orig:
        r = d.parent / (d.name + "_mc2")
        main.append(r if (r / "DONE").exists() else d)
    main = [d for d in main if json.loads((d / "config.json").read_text())["alpha"] == 0.1]
    print("main-study runs (alpha 0.1, _mc2 replaces):", len(main))
    kept_stats(main, "OLD MAIN (atol 1e-5, FD stray field)")
    rule_counts(main, "OLD MAIN")
    a001 = sorted([d for d in RUNS.iterdir() if re.match(r"^C\d_a001$", d.name)])
    kept_stats(a001, "ALPHA 0.01 C#_a001")
    rule_counts(a001, "ALPHA 0.01")
    l6 = sorted([d for d in RUNS.iterdir() if re.match(r"^K1_\d_L6$", d.name)])
    l5 = sorted([d for d in RUNS.iterdir() if re.match(r"^K1_\d$", d.name)])
    kept_stats(l6, "K1_k_L6 (atol 1e-6, main protocol, dx 3)")
    kept_stats(l5, "K1_k (atol 1e-5, same lines)")
    rule_counts(l6, "K1_k_L6"); rule_counts(l5, "K1_k")
    if "--design" in sys.argv:
        compare_rules(main, "OLD MAIN")
