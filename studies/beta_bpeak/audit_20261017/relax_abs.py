import json, sys, numpy as np
sys.path.insert(0, ".")
from scan2 import main_dirs
rows = []
for d in main_dirs():
    s = json.loads((d / "summary.json").read_text()); h = s["init"]["relax_history"]
    E = np.array([x["E_over_Kd"] for x in h]); dE = np.abs(np.diff(E))
    st = [x for x in s["stages"] if x.get("measure")]
    c = st[0]["cycles"][2:] or st[0]["cycles"][-1:]
    w9 = np.mean([x["w_loop"] for x in c if x["w_loop"] > 0])
    rows.append((d.name, abs(E[-1]), dE[-2:].max(), w9, h[-1]["torque_rms"]))
a = np.array([r[1:] for r in rows])
q = lambda v: "min %.3g  10%% %.3g  median %.3g  90%% %.3g  max %.3g" % (v.min(), *np.quantile(v, [.1, .5, .9]), v.max())
print("n", len(a))
print("|E|/(K_d V):", q(a[:, 0]))
print("last |dE| per chunk /(K_d V):", q(a[:, 1]))
print("w_loop(9 mT):", q(a[:, 2]))
print("1e-6 |E| (old relative bound) :", q(1e-6 * a[:, 0]))
for dabs in (1e-10, 3e-10, 1e-9):
    r = 6.67 * dabs / a[:, 2]
    print("delta_abs %.0e: drift per cycle / w9: %s ; runs whose last chunks exceed it: %d" % (dabs, q(r), (a[:, 1] > dabs).sum()))
r = 6.67e-6 * a[:, 0] / a[:, 2]
print("relative 1e-6|E|: drift per cycle / w9:", q(r))
