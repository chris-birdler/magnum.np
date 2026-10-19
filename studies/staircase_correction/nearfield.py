"""
Near-field correction of the stray field for the staircase surface of voxel spheres, after Donahue & McMichael,
IEEE Trans. Magn. 43, 2878 (2007) (PROTOKOLL 3.1 ... 3.4, PLAN.md section 2).

V1 cells (cut cells): a fine mesh with q sub-cells per axis and coarse cell carries the fine body (fine cells whose
centre is inside a sphere). A coarse cell with at least one fine cell of the body is magnetic, with the volume
fraction v_i = n_i / q^3. The coarse FFT field uses the moments v_j m_j; the correction replaces, inside the window
W(i) (max |offset| <= p, minimum image), the coarse coupling v_j K^c(i - j) by the coupling L_ij of the fine cells:

    h_i = h_i(FFT, v_j m_j) + sum_{j in W(i)} (L_ij - v_j K^c(i - j)) m_j            (pairwise form, eq. 5)
    h_i = h_i(FFT, v_j m_j) + [sum_{j in W(i)} (L_ij - v_j K^c(i - j))] m_i          (cell-local form, eqs. 8, 9)

L_ij = (1 / n_i) sum_{a in i} sum_{b in j} K^f(a - b). Sign convention of magnum.np: H = K M (self term about -1/3).
Units: lengths in cells of the respective mesh (the kernel is scale-free), fields h = H / M_s.
"""
import numpy as np
import torch

from magnumnp.field_terms.demag_staircase import (NEWELL_P, SYM, kernel_table, sym_matvec, blocks,  # noqa: F401
                                                  NearCorrection, staircase_cells as v1_cells)


class FFTField:
    """Stray field h = H / M_s of DemagField on one mesh, with the kernel set up once.
    setting: "open" (no PBC), or ("newell", k) / ("newell_ext", k) with k image shells and true_pbc (PROTOKOLL 7.105)."""

    def __init__(self, N, setting, p=NEWELL_P, cache_dir=None):
        from magnumnp import Mesh, State, DemagField
        if setting == "open":
            pbc, term = (0, 0, 0), DemagField(p=p, cache_dir=cache_dir)
        else:
            kind, k = setting
            pbc, term = (k, k, k), DemagField(p=p, true_pbc=True, extrapolate=(kind == "newell_ext"),
                                              cache_dir=cache_dir)
        self.N = N
        self.state = State(Mesh((N, N, N), (1e-9, 1e-9, 1e-9), pbc=pbc))
        self.term = term

    def h(self, ms, m):
        """ms [N, N, N] (moment weight per cell), m [N, N, N, 3] -> h [N, N, N, 3] (double)."""
        st = self.state
        st.material["Ms"] = torch.as_tensor(ms, dtype=st.dtype).reshape(self.N, self.N, self.N, 1)
        st.m = torch.as_tensor(m, dtype=st.dtype)
        with torch.no_grad():
            return self.term.h(st).to(torch.float64)


def fine_owner_m(m_coarse, q):
    """Fine magnetisation of the V1 cells: each fine cell takes m of its coarse cell (the caller masks the body)."""
    m = torch.as_tensor(m_coarse)
    return m.repeat_interleave(q, 0).repeat_interleave(q, 1).repeat_interleave(q, 2)


def coarse_mean(h_fine, occ_fine, q):
    """Mean of the fine field over the body fine cells of each coarse cell (0 where a cell has none)."""
    Nf = h_fine.shape[0]
    N = Nf // q
    w = torch.as_tensor(occ_fine, dtype=torch.float64).reshape(N, q, N, q, N, q)
    hf = h_fine.reshape(N, q, N, q, N, q, 3)
    s = (hf * w[..., None]).sum(dim=(1, 3, 5))
    n = w.sum(dim=(1, 3, 5))
    return s / n.clamp(min=1.0)[..., None]
