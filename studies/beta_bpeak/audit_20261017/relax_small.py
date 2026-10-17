import json, sys, numpy as np
sys.path.insert(0, ".")
from scan2 import main_dirs
rows = []
for d in main_dirs():
    s = json.loads((d / "summary.json").read_text()); h = s["init"]["relax_history"]
    E = np.array([x["E_over_Kd"] for x in h]); dE = np.diff(E)
    rows.append((d.name, abs(E[-1]), s["init"]["relax_iterations"], [h[-2]["rel_dE"], h[-1]["rel_dE"]], np.abs(dE[-2:])))
rows.sort(key=lambda r: r[1])
print("n |E|<1e-5: %d, <1e-4: %d, <3e-4: %d of %d" % tuple([sum(r[1] < x for r in rows) for x in (1e-5, 1e-4, 3e-4)] + [len(rows)]))
print("runs with an exactly zero dE in the last 2 chunks: %d" % sum((r[4] == 0).any() for r in rows))
for r in rows[:8]: print(r[0], "|E| %.2e it %d rel %s |dE| %s" % (r[1], r[2], ["%.1e" % x for x in r[3]], ["%.1e" % x for x in r[4]]))
for k in ("K1_0_mc2","K1_1","K1_2","K2_0","K2_1","K2_2","K3_0","K3_1","K3_2","K4_0_mc2","K4_1_mc2","K4_2"):
    r = [x for x in rows if x[0] == k][0]; print(k, "|E| %.1e it %d" % (r[1], r[2]), end="; ")
print()
