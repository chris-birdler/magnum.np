"""
Geometry of the periodic powder model.

The simulation box is ONE cubic FCC unit cell with edge length a and
periodic boundary conditions in x, y and z. The cell holds 4 spheres
(diameter d) at the FCC lattice sites

    (0, 0, 0), (0, a/2, a/2), (a/2, 0, a/2), (a/2, a/2, 0).

Packing fraction:  phi = 4 * (pi/6) * d^3 / a^3 = (2 pi / 3) * (d/a)^3
Touching spheres:  d = a / sqrt(2)  ->  phi_max = 0.7405

Rules that this module enforces:
  * spheres must not overlap (d < a/sqrt(2)),
  * no two voxels of DIFFERENT particles can touch, not even at an edge or a
    corner (26-neighbourhood): at least one void cell lies between any two
    particles. The exchange stencil of magnum.np couples only face
    neighbours, but an edge or corner contact is not a physical gap and
    gives extreme stray fields.
"""
import math
import numpy as np

__all__ = ["fcc_box", "phi_from", "a_from", "FCC_SITES"]

FCC_SITES = np.array([[0.0, 0.0, 0.0],
                      [0.0, 0.5, 0.5],
                      [0.5, 0.0, 0.5],
                      [0.5, 0.5, 0.0]])


def phi_from(d, a):
    """Packing fraction of 4 spheres with diameter d in a cube with edge a."""
    return 2.0 * math.pi / 3.0 * (d / a) ** 3


def a_from(d, phi):
    """Edge length of the FCC cell for diameter d and packing fraction phi."""
    return d * (2.0 * math.pi / (3.0 * phi)) ** (1.0 / 3.0)


def fcc_box(N, d, dx, sub=0):
    """
    Voxelize 4 spheres on the FCC sites of a periodic cube.

    :param N:  number of cells per edge (the same in x, y, z)
    :param d:  sphere diameter [m]
    :param dx: cell size [m]  (edge length a = N * dx)
    :param sub: 0 = staircase surface (a cell is in the sphere if its centre is);
                > 0 = volume fraction of each surface cell from sub^3 sample points
                ("frac"; cells with frac > 0 belong to the particle)
    :returns: dict with
        "ids"       int8 array [N,N,N], -1 = void, 0..3 = particle index
        "a"         edge length [m]
        "phi_nom"   nominal packing fraction (exact spheres)
        "phi_vox"   packing fraction of the voxel model
        "gap_nom"   nominal surface gap between nearest neighbours [m]
        "d_cells"   sphere diameter in cells
        "frac"      float32 array [N,N,N]: magnetic volume fraction per cell
                    (0 or 1 for sub = 0)
    """
    a = N * dx
    gap = a / math.sqrt(2.0) - d
    if gap <= 0.0:
        raise ValueError("Spheres overlap: d = %.4g m >= a/sqrt(2) = %.4g m" % (d, a / math.sqrt(2.0)))

    # cell centres; the lattice sites sit on cell vertices if N is even
    x = (np.arange(N) + 0.5) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

    ids = -np.ones((N, N, N), dtype=np.int8)
    frac = np.zeros((N, N, N), dtype=np.float32)
    R = 0.5 * d
    r2 = R ** 2
    if sub > 0:
        q = (np.arange(sub) + 0.5) / sub - 0.5                     # sample offsets in a cell (units of dx)
        QX, QY, QZ = [g.reshape(-1) * dx for g in np.meshgrid(q, q, q, indexing="ij")]
    for p, s in enumerate(FCC_SITES * a):
        # minimum-image distance (periodic box)
        ddx = X - s[0]; ddx -= a * np.round(ddx / a)
        ddy = Y - s[1]; ddy -= a * np.round(ddy / a)
        ddz = Z - s[2]; ddz -= a * np.round(ddz / a)
        rr = ddx**2 + ddy**2 + ddz**2
        if sub > 0:
            r = np.sqrt(rr)
            half = 0.5 * math.sqrt(3.0) * dx
            full = r <= R - half
            edge = np.abs(r - R) < half
            ex, ey, ez = ddx[edge][:, None] + QX, ddy[edge][:, None] + QY, ddz[edge][:, None] + QZ
            fe = np.mean(ex**2 + ey**2 + ez**2 <= r2, axis=1).astype(np.float32)
            fp = np.zeros((N, N, N), dtype=np.float32)
            fp[full] = 1.0
            fp[edge] = fe
            inside = fp > 0.0
            if np.any(ids[inside] >= 0):
                raise ValueError("Surface cells of two spheres overlap (particle %d)." % p)
            frac[inside] = fp[inside]
        else:
            inside = rr <= r2
            if np.any(ids[inside] >= 0):
                raise ValueError("Voxelized spheres overlap (particle %d)." % p)
            frac[inside] = 1.0
        ids[inside] = p

    _check_no_contact(ids)

    return {"ids": ids,
            "a": a,
            "phi_nom": phi_from(d, a),
            "phi_vox": float(np.mean(frac)),
            "gap_nom": gap,
            "d_cells": d / dx,
            "frac": frac}


def _check_no_contact(ids):
    """Fail if voxels of different particles touch at a face, an edge or a
    corner (26-neighbourhood, periodic)."""
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for sz in (-1, 0, 1):
                if (sx, sy, sz) == (0, 0, 0):
                    continue
                nb = np.roll(ids, (sx, sy, sz), axis=(0, 1, 2))
                contact = (ids >= 0) & (nb >= 0) & (ids != nb)
                if np.any(contact):
                    raise ValueError("Particles touch on the voxel grid (shift %s, %d cells). "
                                     "Use a smaller d or a finer dx." % ((sx, sy, sz), int(contact.sum())))

