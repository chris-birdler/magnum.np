# box 11: largest |m_r - m_0| at the window ends of a K1 file with repetitions (step_dump --k1_reps, field dm_end), for each window and repetition above 1e-9
import sys
import torch
d = torch.load(sys.argv[1], map_location="cpu", weights_only=False)
n = 0
for k0, w in d["win"].items():
    for r, x in enumerate(w.get("reps", []), 1):
        n += 1
        if x["dm_end"] > 1e-9:
            print("k0 %3d rep %2d: largest |m_r - m_0| %.3e" % (k0, r, x["dm_end"]))
print("%d window ends compared; the others are below 1e-9; largest of them %.3e" % (
    n, max([x["dm_end"] for w in d["win"].values() for x in w.get("reps", []) if x["dm_end"] <= 1e-9], default=0.0)))
