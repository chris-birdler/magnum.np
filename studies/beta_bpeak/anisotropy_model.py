"""
Anisotropy model of the study (fictitious nanocrystalline powder, reduced units).

Two sources, both uniaxial. There is no other anisotropy.

1. Herzer residual anisotropy (grain anisotropy after exchange averaging)
   The grains (nm scale) are not resolved. Herzer: K_eff = K1^4 D^6 / A^3 and
   L_eff = sqrt(A / K_eff). The exchange-averaged anisotropy has a random
   direction on the scale L_eff. Model: cubes with edge L_eff, one random
   isotropic easy axis per cube, amplitude K_eff. In reduced units

       Q_eff = K_eff / K_d = (l_ex / L_eff)^2

   (the cube edge = L_eff is a convention; Herzer fixes L_eff only up to a
   factor of order 1). Translation to a material: Q_eff = Q1 * (D/delta1)^6
   with Q1 = K1/K_d and delta1 = l_ex/sqrt(Q1).

2. Residual stress on the particle scale
   One uniaxial anisotropy K_p = r_p * K_eff per particle. The axes are
   deterministic: cos(theta) = 1/8, 3/8, 5/8, 7/8 to the drive (midpoint rule
   for an isotropic powder), azimuths 0, 90, 180, 270 deg. No orientation
   statistics.

Mesh independence
   The cubes are defined in physical coordinates relative to the particle
   centre (minimum image), with a random offset per particle. The axis table
   depends only on (seed, particle, L_eff, d), not on the mesh. If the box
   edge a, the cube edge L_eff and the offsets are multiples of every cell
   size used (A_QUANT = 12 l_ex: dx = 3, 2 and 1.5 l_ex, particle centres at 0
   and a/2), the cube boundaries lie on cell faces of every mesh, and all
   meshes see exactly the same cubes.
"""
import math
import numpy as np

from geometry import FCC_SITES

__all__ = ["A_QUANT", "CUBE_QUANT", "box_edge", "cube_axes", "particle_axes"]

A_QUANT = 12.0      # l_ex; box edge quantum (particle centres at a/2 -> multiples of 6)
CUBE_QUANT = 6.0    # l_ex; cube edge and offset quantum (common multiple of dx = 3, 2 and 1.5)


def box_edge(d_lex, phi, quant=A_QUANT):
    """
    Box edge a (l_ex) quantized to a multiple of `quant`, and the particle
    diameter that keeps the packing fraction phi exact.
    :returns: (a_lex, d_lex_used)
    """
    a = d_lex * (2.0 * math.pi / (3.0 * phi)) ** (1.0 / 3.0)
    if quant:
        a = quant * max(1, round(a / quant))
    d = a * (3.0 * phi / (2.0 * math.pi)) ** (1.0 / 3.0)
    return a, d


def _cell_centres(N, dx):
    x = (np.arange(N) + 0.5) * dx
    return np.meshgrid(x, x, x, indexing="ij")


def cube_axes(ids, dx, a, d, L, seed):
    """
    Random easy axis per cube of edge L, per particle (all lengths in l_ex).

    :param ids:  int array [N,N,N], -1 = void, 0..3 = particle index (geometry.fcc_box)
    :returns: (axes [N,N,N,3] float64, zero in the void;
               number of cubes that touch each particle [4])
    """
    N = ids.shape[0]
    X, Y, Z = _cell_centres(N, dx)
    axes = np.zeros(ids.shape + (3,))
    n_cubes = []
    M = int(math.ceil((0.5 * d + L) / L)) + 1                 # index range -M..M covers the particle
    for p, s in enumerate(FCC_SITES * a):
        rng = np.random.default_rng([seed, 1000 + p])
        tab = rng.standard_normal((2 * M + 1,) * 3 + (3,))
        tab /= np.linalg.norm(tab, axis=-1, keepdims=True)   # isotropic unit vectors
        off = CUBE_QUANT * rng.integers(0, max(1, int(round(L / CUBE_QUANT))), size=3)
        sel = ids == p
        idx = []
        for c, s_c, o in zip((X, Y, Z), s, off):
            r = c[sel] - s_c
            r -= a * np.round(r / a)                          # minimum image
            idx.append(np.floor((r + o) / L).astype(np.int64) + M)
        axes[sel] = tab[idx[0], idx[1], idx[2]]
        n_cubes.append(len(set(zip(*idx))))
    return axes, n_cubes


def particle_axes(direction):
    """Deterministic easy axes of the 4 particles: cos(theta) = 1/8, 3/8, 5/8, 7/8
    to the drive direction, azimuths 0, 90, 180, 270 deg around the drive."""
    e = np.asarray(direction, float); e /= np.linalg.norm(e)
    t = np.array([0.0, 0.0, 1.0]) if abs(e[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(e, t); e1 /= np.linalg.norm(e1)
    e2 = np.cross(e, e1)
    ax = []
    for c, psi in zip((1/8, 3/8, 5/8, 7/8), (0.0, 0.5 * math.pi, math.pi, 1.5 * math.pi)):
        s = math.sqrt(1.0 - c * c)
        ax.append(c * e + s * (math.cos(psi) * e1 + math.sin(psi) * e2))
    return np.array(ax)
