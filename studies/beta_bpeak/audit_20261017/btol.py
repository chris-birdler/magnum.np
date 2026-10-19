import json, sys, numpy as np
sys.path.insert(0, ".")
from scan2 import main_dirs
dev, devm = [], []
for d in main_dirs():
    s = json.loads((d / "summary.json").read_text())
    for st in s["stages"]:
        if not st.get("measure") or len(st["cycles"]) < 3: continue
        bt = st.get("b_target"); 
        if not bt: continue
        dev.append(st["cycles"][1]["b_peak"] / bt - 1)
        devm.append(np.mean([c["b_peak"] for c in st["cycles"][2:]]) / bt - 1)
dev, devm = np.abs(dev), np.abs(devm)
print("stages", len(dev), "|b(c1)/b_t-1| > 0.01: %.2f, > 0.03: %.2f; median %.4f 90%% %.4f max %.3f" % ((dev > .01).mean(), (dev > .03).mean(), np.median(dev), np.quantile(dev, .9), dev.max()))
print("mean of measured cycles: median %.4f 90%% %.4f max %.3f" % (np.median(devm), np.quantile(devm, .9), devm.max()))
