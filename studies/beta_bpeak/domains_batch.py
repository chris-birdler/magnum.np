#!/usr/bin/env python
"""
Wall measures of all snapshots of one run (domains.py), written to <run>/domains.json.

The run must have m_<stage>_c<cycle>_ph<deg>.vti (--snap_phases 4, or more: then the phases 90, 180, 270 and
360 deg are taken; geometry.vti is used if it exists). Per stage: the 4 phases of the LAST snapshot cycle
(H max, H = 0 falling, H min, H = 0 rising), B_peak of the stage from summary.json. init.pt (the relaxed
virgin state, key "init") and checkpoint.pt (the state at the end of the run, key "final") are analysed
too if they exist.

    python domains_batch.py runs/B300_s921 [runs/...]
"""
import json
import pathlib
import re
import sys

import numpy as np

import domains

JS_T = 1.5


def select(d):
    """Snapshot files per stage: the last snapshot cycle only (never mix cycles); of it the phases 90, 180,
    270, 360 deg if all exist, else all phases of that cycle in order."""
    groups = {}
    for f in sorted(pathlib.Path(d).glob("m_*_c*_ph*.vti")):
        m = re.match(r"m_(.+)_c(\d+)_ph(\d+)\.vti", f.name)
        groups.setdefault(m.group(1), {}).setdefault(int(m.group(2)), {})[int(m.group(3))] = f
    sel = {}
    for name, cyc in groups.items():
        ph = cyc[max(cyc)]
        quarter = (90, 180, 270, 360)
        sel[name] = [ph[k] for k in quarter] if all(k in ph for k in quarter) else [ph[k] for k in sorted(ph)]
    return sel


def run(d):
    d = pathlib.Path(d)
    summ = json.loads((d / "summary.json").read_text())
    b_of = {}
    for st in summ["stages"]:
        if st.get("cycles"):
            cy = st["cycles"][2:] or st["cycles"][-1:]
            b_of[st["name"]] = float(np.mean([c["b_peak"] for c in cy])) * JS_T * 1000.0
    out = {"run": d.name, "stages": []}
    for key, fname in (("init", "init.pt"), ("final", "checkpoint.pt")):
        if (d / fname).exists():
            out[key] = domains.measures(domains.load_state(d / fname))
    for name, fs in select(d).items():
        res = domains.analyze([domains.load_state(f) for f in fs])
        res.update({"stage": name, "B_peak_mT": b_of.get(name, float("nan")),
                    "files": [f.name for f in fs]})
        out["stages"].append(res)
    out["stages"].sort(key=lambda r: r["B_peak_mT"])
    (d / "domains.json").write_text(json.dumps(out, indent=1))
    print("%s: %d stages%s%s" % (d.name, len(out["stages"]), ", init" if "init" in out else "",
                                 ", final" if "final" in out else ""), flush=True)


if __name__ == "__main__":
    for d in sys.argv[1:]:
        run(d)
