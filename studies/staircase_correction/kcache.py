"""
Kernel cache of DemagField for the tests of PLAN 14: run_loops.py gives DemagField no cache folder, so each set-up
computes the Newell kernel again. use(folder) gives each DemagField without a folder this one. The cache name has the
mesh and the extrapolation, but no dtype, and demag.py casts a cached kernel of another dtype; thus use one folder for
each geometry and precision. use() stops the run if a cached kernel has another dtype than the run.
"""
import os

import torch


def use(folder):
    from magnumnp.field_terms.demag import DemagField
    f_init, f_init_N = DemagField.__init__, DemagField._init_N

    def init(self, *a, **k):
        f_init(self, *a, **k)
        if self._cache_dir is None:
            self._cache_dir = folder

    def init_N(self, state):
        name = "/N_%s%s.pt" % (str(state.mesh).replace(" ", ""), "_ext" if self._extrapolate else "")
        p = self._cache_dir + name if self._cache_dir is not None else None
        if p is not None and os.path.isfile(p):
            dt = torch.load(p, map_location="cpu")[0].dtype
            if dt != torch.get_default_dtype():
                raise SystemExit("kcache: %s has dtype %s, the run %s (one folder for each precision)"
                                 % (p, dt, torch.get_default_dtype()))
        return f_init_N(self, state)

    DemagField.__init__, DemagField._init_N = init, init_N
