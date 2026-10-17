import json, sys, numpy as np
sys.path.insert(0, ".")
from scan2 import main_dirs, RUNS
rel, cyc, tot, ncyc, names = [], [], [], [], []
for d in main_dirs():
    s = json.loads((d / "summary.json").read_text())
    r = s["init"].get("wall_s", 0.0)
    cw = [c.get("wall_s", np.nan) for st in s["stages"] for c in st["cycles"]]
    rel.append(r); cyc.append(np.nansum(cw)); ncyc.append(len(cw)); names.append(d.name)
rel, cyc, ncyc = map(np.array, (rel, cyc, ncyc))
print("runs", len(rel), "cycles per run median", np.median(ncyc))
print("relax h: median %.2f  sum %.1f ; cycles h: median %.2f  sum %.1f ; relax share %.2f" % (
    np.median(rel)/3600, rel.sum()/3600, np.median(cyc)/3600, cyc.sum()/3600, rel.sum()/(rel.sum()+cyc.sum())))
for k in ("K1_0","K2_0","K3_0","K4_0"):
    if k in names:
        i = names.index(k); print(k, "relax h %.2f cycles h %.2f n %d" % (rel[i]/3600, cyc[i]/3600, ncyc[i]))
k = [i for i, n in enumerate(names) if n.startswith("K")]
print("corners", [names[i] for i in k])
print("corner relax h sum %.1f, cycles h sum %.1f, max run h %.2f" % (rel[k].sum()/3600, cyc[k].sum()/3600, max((rel[i]+cyc[i])/3600 for i in k)))
