import sys, numpy as np
sys.path.insert(0, "/home/christoph/Documents/code/magnum.np-newell/studies/beta_bpeak")
sys.path.insert(0, ".")
import eval_corners_nx as E, analyze_design as AD
from scan2 import RUNS
for k in "1234":
    b = []
    for s in "012":
        d = RUNS / ("K%s_%s_mc2" % (k, s)); d = d if (d / "DONE").exists() else RUNS / ("K%s_%s" % (k, s))
        y = E.lnw(d); x = np.log(np.array(AD.AMPS_MT)); m = np.isfinite(y)
        b.append(np.polyfit(x[m], y[m], 1)[0])
    print("K%s beta per seed %s  SD %.3f" % (k, np.round(b, 3), np.std(b, ddof=1)))
