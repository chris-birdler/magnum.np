import json, sys, numpy as np
sys.path.insert(0, ".")
from scan2 import main_dirs
sys.path.insert(0, "/home/christoph/Documents/code/magnum.np-newell/studies/beta_bpeak")
import run_loops as RL
tot, extra, short, nst, nmc2, maxc, miss = 0, 0, 0, 0, 0, 0, 0
for d in main_dirs():
    s = json.loads((d / "summary.json").read_text())
    for st in s["stages"]:
        if not st.get("measure"): continue
        n = len(st["cycles"]); tot += n; nst += 1
        if d.name.endswith("_mc2"):
            extra += max(0, n - 7); maxc = max(maxc, n)
        ct = RL.calm_tail(st["cycles"][2:7], 0.20)
        short += ct < 2
        miss += RL.calm_tail(st["cycles"][2:], 0.20) < 2
    nmc2 += d.name.endswith("_mc2")
print("stages %d, cycles %d; stages with calm tail < 2 after 7 cycles: %d (%.1f %%); _mc2 runs %d, extra cycles %d (%.1f %% of all cycles), max cycles of a stage %d; stages still < 2 calm at the end: %d" % (nst, tot, short, 100*short/nst, nmc2, extra, 100*extra/tot, maxc, miss))
