#!/usr/bin/env python
"""
Time and memory of the V1 tie (thin_face, D0 of PLAN 4.21) at K3 of jobs/plan421.txt at q 16, in one process with one
thread, written to results/plan421_tie_time.json: the wall time (s) of geometry.fcc_box (the fraction surface), of
geometry.fcc_subcells (the sub-cells, as run_loops) and of field_terms_extra.tie_thin_face, the tie counts, and the
peak resident memory of the process (GiB). The meshes of K2 and K4 (152^3 and 156^3) are not measured.
    OMP_NUM_THREADS=1 python time_tie421.py
"""
import json
import pathlib
import resource
import shlex
import time

import scan421
from field_terms_extra import tie_thin_face
from geometry import fcc_box, fcc_subcells

HERE = pathlib.Path(__file__).resolve().parent

if __name__ == "__main__":
    line = [l for l in (HERE / "jobs" / "plan421.txt").read_text().splitlines() if l.startswith("K3_0_v16 ")][0]
    N, d, dx, _, _ = scan421.mesh(line.split(None, 1)[1])
    args = shlex.split(line.split(None, 1)[1])
    dx_lex = float(args[args.index("--dx_lex") + 1])
    res = {"run": "K3_0_v16", "N": int(N), "q": 16}
    t = time.time()
    b = fcc_box(N, d, dx, sub=16)
    res["fcc_box_s"] = time.time() - t
    t = time.time()
    sub = fcc_subcells(N, d, dx, 16)[2]
    res["fcc_subcells_s"] = time.time() - t
    t = time.time()
    tm = tie_thin_face(b["ids"], b["frac"], dx_lex, sub, 16)
    res["tie_thin_face_s"] = time.time() - t
    res["tie"] = tm.info
    res["peak_rss_GiB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
    (HERE / "results" / "plan421_tie_time.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
