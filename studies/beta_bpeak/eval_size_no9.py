#!/usr/bin/env python
"""
Mesh check without the 9 mT point (request of Chris 2026-10-12, after the M212 result; post hoc, not pre-registered):
the 9 mT point is suspect (the core moves less than one coarse cell, 7.74). For the three measured points (K3, base
microstructure at d/l_ex 212 and 300): β from 50 and 100 mT only (slope of ln P over ln B per seed), the mesh error
Δβ = β(dx 3) − β(dx 1.5) paired per seed, and the size lever 212 → 300 on each mesh (Welch, as eval_size.py). For K3
also β from 20 … 100 mT (5 amplitudes, eval_k3fill runs). Same statistics as eval_size.py (SE, 90 % CI with t).

    python eval_size_no9.py > results/mesh_no9.md
"""
import numpy as np

import eval_k3fill as KF
import eval_size as ES

AMPS2 = (50, 100)


def betas(group, amps):
    """{mesh: {seed: β}} from the stage dicts of eval_size.load_group."""
    out = {"dx15": {}, "dx3": {}}
    for (mesh, s), st in group.items():
        b = ES.beta(st, amps)
        if b is not None:
            out[mesh][s] = b
    return out


def k3_betas(amps):
    """β of the K3 mesh pairs (eval_k3fill folders) over the given amplitudes, {mesh: {seed: β}}."""
    out = {"dx15": {}, "dx3": {}}
    for mesh, dx in (("dx15", "15"), ("dx3", "3")):
        for s in KF.SEEDS:
            st = {}
            for b in amps:                    # one folder per amplitude: its single steady stage (eval_size.nominal
                v = list(ES.stages(KF.folder(dx, s, b)).values())       # knows only 9 / 50 / 100 mT)
                if len(v) == 1:
                    st[b] = v[0]
            if all(b in st for b in amps):
                out[mesh][s] = np.polyfit(np.log([st[b][1] for b in amps]), np.log([st[b][0] for b in amps]), 1)[0]
    return out


def main():
    L = ["# Mesh check without the 9 mT point (post hoc, request of Chris 2026-10-12)", "",
         "β from 50 and 100 mT only (and for K3 also from 20 … 100 mT). dx in l_ex. ± SE, 90 %% CI with t. Fewer than "
         "%d seeds: no statement. The full table with 9 mT: results/mesh212.md." % ES.MIN_N, ""]
    groups = {}
    for prefix, lab, seeds in ES.GROUPS:
        groups[prefix] = (lab, betas(ES.load_group(prefix, seeds), AMPS2))
    L += ["## β (50 → 100 mT) per mesh and mesh error", "",
          "| point | β fine dx 1.5 | β coarse dx 3 | Δβ = coarse − fine (paired) |", "|---|---|---|---|"]
    for prefix, (lab, b) in groups.items():
        paired = [b["dx3"][s] - b["dx15"][s] for s in b["dx15"] if s in b["dx3"]]
        L.append("| %s | %s | %s | %s |" % (lab, ES.fmt(ES.summary(b["dx15"].values()), n=len(b["dx15"])),
                                             ES.fmt(ES.summary(b["dx3"].values()), n=len(b["dx3"])),
                                             ES.fmt(ES.summary(paired), n=len(paired))))
    for amps, tag in (((20, 35, 50, 70, 100), "20 … 100 mT"), (AMPS2, "50 → 100 mT, eval_k3fill runs")):
        b = k3_betas(amps)
        paired = [b["dx3"][s] - b["dx15"][s] for s in b["dx15"] if s in b["dx3"]]
        L.append("| K3, %s | %s | %s | %s |" % (tag, ES.fmt(ES.summary(b["dx15"].values()), n=len(b["dx15"])),
                                                 ES.fmt(ES.summary(b["dx3"].values()), n=len(b["dx3"])),
                                                 ES.fmt(ES.summary(paired), n=len(paired))))
    L += ["", "## Size lever 212 → 300 (base microstructure; value at 1 µm minus at 0.71 µm, Welch)", "",
          "| output | lever fine dx 1.5 | lever coarse dx 3 | mesh share (coarse − fine) |", "|---|---|---|---|"]
    b212, b300 = groups["M212"][1], groups["M300"][1]
    lev = {}
    for mesh in ("dx15", "dx3"):
        lev[mesh] = ES.difference(ES.summary(b300[mesh].values()), ES.summary(b212[mesh].values()))
    f = lambda d: "–" if d is None else "%+.2f ± %.2f (CI %+.2f … %+.2f)" % (d[0], d[1], d[0] - d[2], d[0] + d[2])
    e300 = ES.summary([b300["dx3"][s] - b300["dx15"][s] for s in b300["dx15"] if s in b300["dx3"]])
    e212 = ES.summary([b212["dx3"][s] - b212["dx15"][s] for s in b212["dx15"] if s in b212["dx3"]])
    L.append("| β (50 → 100 mT) | %s | %s | %s |" % (f(lev["dx15"]), f(lev["dx3"]), f(ES.difference(e300, e212))))
    print("\n".join(L))
    import json
    import pathlib
    js = {"points": {}, "lever": {}}                        # results/mesh_no9.json for the deck (added 2026-10-13)
    for prefix, (lab, b) in groups.items():
        paired = [b["dx3"][s] - b["dx15"][s] for s in b["dx15"] if s in b["dx3"]]
        sm = ES.summary(paired)
        js["points"][prefix] = {"label": lab, "dbeta": ({"mean": sm[0], "se": sm[1], "n": sm[2], "lo": sm[0] - sm[3],
                                                          "hi": sm[0] + sm[3]} if sm else None),
                                "beta": {m: list(b[m].values()) for m in ("dx15", "dx3")}}
    for m, d in (("dx15", lev["dx15"]), ("dx3", lev["dx3"]), ("share", ES.difference(e300, e212))):
        js["lever"][m] = {"mean": d[0], "se": d[1], "lo": d[0] - d[2], "hi": d[0] + d[2]} if d else None
    (pathlib.Path(__file__).resolve().parent / "results" / "mesh_no9.json").write_text(json.dumps(js, indent=1))


if __name__ == "__main__":
    main()
