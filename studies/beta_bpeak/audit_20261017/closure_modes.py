import json, re, sys, numpy as np
sys.path.insert(0, ".")
from scan2 import main_dirs, RUNS
from scan_rules import stage_cycles, make_rule, eps
base = make_rule()

def stage_M(d, si, name):
    f = d / ("samples_%02d_%s.csv" % (si, name))
    a = np.loadtxt(f, delimiter=",", comments="#", usecols=(0, 3), ndmin=2)
    out = {}
    for c in np.unique(a[:, 0]).astype(int):
        m = a[a[:, 0] == c, 1]
        out[c] = (m[0], m[-1])
    return out

def analyse(dirs, label):
    res = []
    for d in dirs:
        s = json.loads((d / "summary.json").read_text())
        k = 0
        for si, st in enumerate(s["stages"]):
            if not st.get("measure") or not st.get("cycles"):
                continue
            fm = 3 if "h_amp_second" in st else (2 if "h_amp_initial" in st else 1)
            t = base(st["cycles"][fm:], None)
            if len(t) < 3:
                k += 1; continue
            M = stage_M(d, si, st["name"])
            Mpk = np.mean([c["M_peak"] for c in t])
            starts = [M[c["cycle"]][0] for c in t] + [M[t[-1]["cycle"]][1]]   # M at H = 0 before each kept cycle and after the last
            x = np.array(starts) / Mpk
            dM = np.diff(x)                       # signed closure per cycle
            r2 = np.array([abs(x[i + 2] - x[i]) / max(0.5 * (abs(dM[i]) + abs(dM[i + 1])), 1e-30) for i in range(len(x) - 2)])
            same = np.mean(np.sign(dM[1:]) == np.sign(dM[:-1]))
            res.append(dict(run=d.name, amp=k, clo=np.mean(np.abs(dM)), r2=np.median(r2), same=same, net=abs(x[-1] - x[0]) / (len(dM) * np.mean(np.abs(dM))),
                            em=sum(c["w_dis"] for c in t) / sum(c["w_loop"] for c in t) - 1, n=len(t)))
            k += 1
    r2 = np.array([r["r2"] for r in res]); clo = np.array([r["clo"] for r in res]); net = np.array([r["net"] for r in res])
    print("\n## %s: %d stages with >= 3 kept cycles" % (label, len(res)))
    print("2-cycle/1-cycle closure ratio r2 (iid scatter ~1, period-2 ~0, drift ~2): quantiles 10/50/90 %s" % np.round(np.quantile(r2, [.1, .5, .9]), 2))
    print("net/total path (drift ~1, back-and-forth ~0): quantiles 10/50/90 %s" % np.round(np.quantile(net, [.1, .5, .9]), 2))
    for lo, hi in ((0, 3e-4), (3e-4, 1e-3), (1e-3, 2e-3), (2e-3, 1)):
        m = (clo >= lo) & (clo < hi)
        if m.sum():
            print("  closure %.0e..%.0e: n %3d  r2<0.5 (period-2-like) %3d  r2>1.6 (drift-like) %3d  net>0.8 %3d  median net %.2f" % (
                lo, hi, m.sum(), (r2[m] < 0.5).sum(), (r2[m] > 1.6).sum(), (net[m] > 0.8).sum(), np.median(net[m])))
    return res

r = analyse(main_dirs(), "OLD MAIN alpha 0.1")
a = analyse(sorted(d for d in RUNS.iterdir() if re.match(r"^C\d_a001$", d.name)), "ALPHA 0.01 (baseline tail >= 3 cycles only)")
# alpha 0.01: all measured cycles, signed eps and signed dM pattern, no tail filter
print("\n## ALPHA 0.01, all measured cycles: signed eps sign balance and lag-1 autocorrelation of eps per stage")
ac, pos, neg = [], 0, 0
for d in sorted(d for d in RUNS.iterdir() if re.match(r"^C\d_a001$", d.name)):
    s, sts = stage_cycles(d)
    for st, cyc in sts:
        e = np.array([eps(c) for c in cyc if c["w_loop"] > 0])
        pos += (e > 0.2).sum(); neg += (e < -0.2).sum()
        if len(e) > 3:
            ac.append(np.corrcoef(e[:-1], e[1:])[0, 1])
print("cycles with eps > +0.2: %d, eps < -0.2: %d; median lag-1 autocorr of eps %.2f (negative: alternating carry-over)" % (pos, neg, np.nanmedian(ac)))
