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
  * no two voxels of DIFFERENT particles can be face neighbours
    (else the exchange term couples the particles through the contact).
"""
import math
import numpy as np

__all__ = ["fcc_box", "phi_from", "a_from", "random_rotations", "FCC_SITES"]

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


def fcc_box(N, d, dx):
    """
    Voxelize 4 spheres on the FCC sites of a periodic cube.

    :param N:  number of cells per edge (the same in x, y, z)
    :param d:  sphere diameter [m]
    :param dx: cell size [m]  (edge length a = N * dx)
    :returns: dict with
        "ids"       int8 array [N,N,N], -1 = void, 0..3 = particle index
        "a"         edge length [m]
        "phi_nom"   nominal packing fraction (exact spheres)
        "phi_vox"   packing fraction of the voxel model
        "gap_nom"   nominal surface gap between nearest neighbours [m]
        "d_cells"   sphere diameter in cells
    """
    a = N * dx
    gap = a / math.sqrt(2.0) - d
    if gap <= 0.0:
        raise ValueError("Spheres overlap: d = %.4g m >= a/sqrt(2) = %.4g m" % (d, a / math.sqrt(2.0)))

    # cell centres; the lattice sites sit on cell vertices if N is even
    x = (np.arange(N) + 0.5) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

    ids = -np.ones((N, N, N), dtype=np.int8)
    r2 = (0.5 * d) ** 2
    for p, s in enumerate(FCC_SITES * a):
        # minimum-image distance (periodic box)
        ddx = X - s[0]; ddx -= a * np.round(ddx / a)
        ddy = Y - s[1]; ddy -= a * np.round(ddy / a)
        ddz = Z - s[2]; ddz -= a * np.round(ddz / a)
        inside = ddx**2 + ddy**2 + ddz**2 <= r2
        if np.any(ids[inside] >= 0):
            raise ValueError("Voxelized spheres overlap (particle %d)." % p)
        ids[inside] = p

    _check_no_contact(ids)

    return {"ids": ids,
            "a": a,
            "phi_nom": phi_from(d, a),
            "phi_vox": float(np.mean(ids >= 0)),
            "gap_nom": gap,
            "d_cells": d / dx}


def _check_no_contact(ids):
    """Fail if voxels of different particles are face neighbours (periodic)."""
    for axis in range(3):
        nb = np.roll(ids, 1, axis=axis)
        contact = (ids >= 0) & (nb >= 0) & (ids != nb)
        if np.any(contact):
            raise ValueError("Particles touch on the voxel grid (axis %d, %d faces). "
                             "Use a smaller d or a finer dx." % (axis, int(contact.sum())))


def random_rotations(n, seed):
    """
    n uniformly distributed random rotation matrices (Haar measure).
    Row vectors R[i,0], R[i,1], R[i,2] are the cubic axes of crystal i.
    """
    rng = np.random.default_rng(seed)
    R = np.empty((n, 3, 3))
    for i in range(n):
        q, r = np.linalg.qr(rng.standard_normal((3, 3)))
        q = q * np.sign(np.diag(r))       # unique QR -> Haar distribution
        if np.linalg.det(q) < 0.0:
            q[:, 0] = -q[:, 0]            # proper rotation
        R[i] = q.T
    return R
