# box 11: find the tries of a sample interval n with a dt near the given dt in K1 files (norm and acceptance)
import sys, torch
n0, dt0 = int(sys.argv[1]), float(sys.argv[2])
for fn in sys.argv[3:]:
    d = torch.load(fn, map_location="cpu", weights_only=False)
    f = d["fields"]; i_n, i_dt, i_no = f.index("n"), f.index("dt"), f.index("norm")
    i_acc = f.index("accept") if "accept" in f else None
    for k, w in d["win"].items():
        reps = w["rows"] if "rows" in w else []
        for j, r in enumerate(reps):
            if r[i_n] == n0 and abs(r[i_dt] - dt0) < 2e-3 * dt0:
                print(fn.split("/")[-1], "win", k, "row", j, "rep", r[f.index("rep")] if "rep" in f else "-", "dt %.6e norm %.6f" % (r[i_dt], r[i_no]), "accept", r[i_acc] if i_acc is not None else "-")
