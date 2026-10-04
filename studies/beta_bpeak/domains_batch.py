#!/usr/bin/env python
"""
Wall measures of all snapshots of one run (domains.py), written to <run>/domains.json.

The run must have m_<stage>_c<cycle>_ph<deg>.vti (--snap_phases 4) and geometry.vti; init.pt
(the relaxed virgin state) is analysed too if it exists. Per stage: the 4 phases of the last
cycle (H max, H = 0 falling, H min, H = 0 rising), B_peak of the stage from summary.json.

    python domains_batch.py runs/B300_s921 [runs/...]
"""
import json
import pathlib
import re
import sys

import numpy as np

import domains

JS_T = 1.5


def run(d):
    d = pathlib.Path(d)
    summ = json.loads((d / "summary.json").read_text())
    b_of = {}
    for st in summ["stages"]:
        if st.get("cycles"):
            cy = st["cycles"][2:] or st["cycles"][-1:]
            b_of[st["name"]] = float(np.mean([c["b_peak"] for c in cy])) * JS_T * 1000.0
    groups = {}
    for f in sorted(d.glob("m_*_c*_ph*.vti")):
        m = re.match(r"m_(.+)_c(\d+)_ph(\d+)\.vti", f.name)
        groups.setdefault(m.group(1), []).append((int(m.group(3)), f))
    out = {"run": d.name, "stages": []}
    if (d / "init.pt").exists():
        out["init"] = domains.measures(domains.load_state(d / "init.pt"))
    for name, fs in groups.items():
        fs = [f for _, f in sorted(fs)]
        res = domains.analyze([domains.load_state(f) for f in fs])
        res.update({"stage": name, "B_peak_mT": b_of.get(name, float("nan")),
                    "files": [f.name for f in fs]})
        out["stages"].append(res)
    out["stages"].sort(key=lambda r: r["B_peak_mT"])
    (d / "domains.json").write_text(json.dumps(out, indent=1))
    print("%s: %d stages%s" % (d.name, len(out["stages"]), ", init" if "init" in out else ""), flush=True)


if __name__ == "__main__":
    for d in sys.argv[1:]:
        run(d)
