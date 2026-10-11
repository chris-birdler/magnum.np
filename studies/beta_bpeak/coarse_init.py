#!/usr/bin/env python
"""
Same-state mesh test (PROTOKOLL 7.63, go Chris 2026-10-11): coarse-grain the relaxed virgin state of a dx 1.5 job to
the dx 3 grid, so that dx 3 can start from the state that dx 1.5 found.
  m_coarse(cell) = normalised mean of m over the 2 × 2 × 2 fine cells that lie inside a particle on the fine grid
  (the dx 3 grid has exactly half the cells per edge and the same box; checked). A coarse particle cell without a
  fine particle cell (surface steps) keeps m of the coarse template.
The output keeps "config" and "init" of the coarse template (the relaxed dx 3 state of the same seed), so that
run_loops --init_from accepts it for the dx 3 job arguments; "init" gets the source of m.

    python coarse_init.py FINE_DIR COARSE_TEMPLATE_DIR OUT_FILE
      FINE_DIR             a dx 1.5 run folder with init.pt and config.json (e.g. runs/MK3_dx15_s2031_b9)
      COARSE_TEMPLATE_DIR  the dx 3 run folder of the same seed with init.pt and config.json (e.g. runs/MK3_dx3_s2031_b9)
"""
import json
import pathlib
import sys

import numpy as np
import torch

from geometry import fcc_box


def ids_of(cfg):
    return fcc_box(int(cfg["N_used"]), float(cfg["d"]), float(cfg["dx"]))["ids"]


def coarse_grain(m_f, ids_f, ids_c, m_template):
    """2 × 2 × 2 mean of m over the fine particle cells, normalised; returns m_c and statistics."""
    Nc = ids_c.shape[0]
    if ids_f.shape[0] != 2 * Nc:
        raise ValueError("fine grid %d is not 2 × coarse grid %d" % (ids_f.shape[0], Nc))
    w = (ids_f >= 0).astype(np.float64)
    s = (m_f * w[..., None]).reshape(Nc, 2, Nc, 2, Nc, 2, 3).sum(axis=(1, 3, 5))
    n = w.reshape(Nc, 2, Nc, 2, Nc, 2).sum(axis=(1, 3, 5))
    norm = np.linalg.norm(s, axis=-1)
    mag_c = ids_c >= 0
    ok = (n > 0) & (norm > 0)
    m_c = np.array(m_template, dtype=np.float64)
    m_c[ok] = s[ok] / norm[ok][:, None]
    stats = {"coarse_particle_cells": int(mag_c.sum()), "kept_from_template": int((mag_c & ~ok).sum()),
             "mean_alignment": float((norm[ok & mag_c] / n[ok & mag_c]).mean())}
    return m_c, stats


def main():
    fine, tmpl, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
    cf, cc = (json.loads((d / "config.json").read_text()) for d in (fine, tmpl))
    if abs(cf["a_lex"] - cc["a_lex"]) > 1e-9 or abs(cf["d"] - cc["d"]) > 1e-15 or cf["seed"] != cc["seed"]:
        raise ValueError("box, particle size or seed differ between %s and %s" % (fine, tmpl))
    F = torch.load(fine / "init.pt", map_location="cpu")
    Tm = torch.load(tmpl / "init.pt", map_location="cpu")
    m_c, st = coarse_grain(F["m"].double().numpy(), ids_of(cf), ids_of(cc), Tm["m"].double().numpy())
    init = dict(Tm["init"])
    init["m_source"] = "2x2x2 mean of %s/init.pt (dx %g), coarse_init.py" % (fine.name, cf["dx_lex"])
    init["coarse_grain"] = st
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"m": torch.tensor(m_c, dtype=Tm["m"].dtype), "init": init, "config": Tm["config"]}, out)
    print("%s: %s" % (out, json.dumps(st)))


if __name__ == "__main__":
    main()
