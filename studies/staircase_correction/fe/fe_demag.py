#!/usr/bin/env python
"""
FE comparison, preparation (PROTOKOLL 3.24): stray-field energy of a coarse state with magnum.pi (FE, P1 potential,
periodic box, mean field zero). Runs in the magnum.pi container:

    docker run --rm -v ~/Documents/code/magnum.pi:/home/firedrake/magnum.pi -v WORK:/work magnum.pi.dev-3.13.3 \
        bash -c 'cd /work && python3 fe_demag.py --data fe_data/X3_s951 --mesh box:48 --body vox'

--data   folder of fe_states.py (m.npy, ids.npy, meta.json)
--mesh   box:n  structured mesh of n^3 cubes of the box, 6 tetrahedra each (conforms to the FD cells if n is a multiple
                of N); FILE.msh  gmsh mesh of the box (fe_mesh.py)
--body   vox     the voxel body of the FD mesh (cell centre in a sphere): FD Newell is exact for it
         q<k>    the staircase of the sub-cells dx/k (X3 of the FD side is exact for it)
         full    the whole box with the m of m_full.npy (no surface; FD exact: E_full)
         smooth  the smooth spheres
--m      state   the coarse m of the state; interp: a continuous m from it; 100, 111: uniform m (target of the uniform case: φ (1 − φ) / 3 K_d V_box)

For a P1 potential the gradient is constant in a tetrahedron, so the load and the energy need only the mean of
M = M_s m over each tetrahedron (DG0). For P2 the gradient is linear, so the L2 projection of M onto P1 in each
tetrahedron (DG1, nodes at the vertices) gives the exact load and energy. Both come from the centroids of the 8^L
sub-tetrahedra of L midpoint refinements (equal volumes; --levels L), so they are exact only up to this sampling where
M is not constant in a tetrahedron. The FE energy of the Galerkin solution is not larger than the exact energy for
the sampled sources.
"""
import argparse
import json
import pathlib
import resource
import time

import numpy as np
import firedrake as fd
from firedrake.petsc import PETSc
import magnumpi as mp
from magnumpi.lib import PeriodicMesh

MU0 = 4e-7 * np.pi


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def ref_points(levels):
    """Barycentric coordinates of the centroids of the 8^levels sub-tetrahedra of midpoint refinement."""
    tets = [np.eye(4)]
    for _ in range(levels):
        new = []
        for t in tets:
            v = list(t)
            mid = {(i, j): 0.5 * (v[i] + v[j]) for i in range(4) for j in range(i + 1, 4)}
            m = lambda i, j: mid[(min(i, j), max(i, j))]  # noqa: E731
            new += [np.array([v[0], m(0, 1), m(0, 2), m(0, 3)]), np.array([m(0, 1), v[1], m(1, 2), m(1, 3)]),
                    np.array([m(0, 2), m(1, 2), v[2], m(2, 3)]), np.array([m(0, 3), m(1, 3), m(2, 3), v[3]]),
                    # the inner octahedron, cut along the diagonal m(0,2)-m(1,3)
                    np.array([m(0, 1), m(0, 2), m(0, 3), m(1, 3)]), np.array([m(0, 1), m(0, 2), m(1, 2), m(1, 3)]),
                    np.array([m(0, 2), m(0, 3), m(1, 3), m(2, 3)]), np.array([m(0, 2), m(1, 2), m(1, 3), m(2, 3)])]
        tets = new
    return np.array([t.mean(0) for t in tets])


class Source:
    """M = M_s m at points x (periodic box [0, a)^3), for the voxel body or the smooth spheres."""
    def __init__(self, data, body, mode):
        meta = json.loads((data / "meta.json").read_text())
        self.N, self.dx, self.a, self.R = meta["N"], meta["dx"], meta["a"], 0.5 * meta["d"]
        self.Ms, self.sites = meta["Ms"], np.array(meta["sites"])
        self.ids = np.load(data / "ids.npy")
        self.mode = mode
        if mode in ("state", "interp"):
            self.m = np.load(data / ("m_full.npy" if body == "full" else "m.npy"))
        else:
            v = np.array([1.0, 0.0, 0.0] if mode == "100" else [1.0, 1.0, 1.0])
            self.m = np.broadcast_to(v / np.linalg.norm(v), (self.N,) * 3 + (3,))
        self.body = body
        self.meta = meta
        if body.startswith("q"):
            self.occ = np.load(data / ("occ_%s.npy" % body))

    def cut_cells(self, verts):
        """Cells that are not inside one sphere with all 4 vertices and do not lie outside all spheres with their
        bounding sphere: the surface can cut them."""
        cen = verts.mean(1)
        rb = np.sqrt(((verts - cen[:, None]) ** 2).sum(-1)).max(1)
        inside = np.zeros(len(verts), dtype=bool)
        near = np.zeros(len(verts), dtype=bool)
        for s in self.sites:
            d = cen - s
            shift = self.a * np.round(d / self.a)
            dv = verts - s - shift[:, None]
            inside |= ((dv * dv).sum(-1) < self.R * self.R).all(1)
            near |= np.sqrt(((d - shift) ** 2).sum(1)) - rb < self.R
        return near & ~inside

    def __call__(self, x):
        x = np.mod(x, self.a)
        c = np.minimum((x / self.dx).astype(np.int64), self.N - 1)
        if self.body == "full":                       # the whole box, m from m_full.npy (no surface)
            inside = np.ones(len(x), dtype=bool)
        elif self.body == "vox":
            inside = self.ids[c[:, 0], c[:, 1], c[:, 2]] >= 0
        elif self.body.startswith("q"):                 # body of the q-staircase (each sub-cell takes m of its cell)
            k = int(self.body[1:])
            f = np.minimum((x / (self.dx / k)).astype(np.int64), k * self.N - 1)
            inside = self.occ[f[:, 0], f[:, 1], f[:, 2]]
        else:
            inside = np.zeros(len(x), dtype=bool)
            for s in self.sites:
                d = x - s
                d -= self.a * np.round(d / self.a)
                inside |= (d * d).sum(1) < self.R * self.R
        if self.mode == "interp":                   # continuous m: trilinear in the cell centres, unit length
            g = x / self.dx - 0.5
            i0 = np.floor(g).astype(np.int64)
            w = g - i0
            acc = np.zeros((len(x), 3))
            wsum = np.zeros(len(x))
            for di in (0, 1):
                for dj in (0, 1):
                    for dk in (0, 1):
                        ii = [(i0[:, 0] + di) % self.N, (i0[:, 1] + dj) % self.N, (i0[:, 2] + dk) % self.N]
                        wt = (np.where(di, w[:, 0], 1 - w[:, 0]) * np.where(dj, w[:, 1], 1 - w[:, 1])
                              * np.where(dk, w[:, 2], 1 - w[:, 2]))
                        mc = self.m[ii[0], ii[1], ii[2]]
                        has = (mc * mc).sum(1) > 0               # only cells with m (normalised weights)
                        acc += (wt * has)[:, None] * mc
                        wsum += wt * has
            mi = acc / np.maximum(wsum, 1e-30)[:, None]
            M = self.Ms * mi / np.maximum(np.linalg.norm(mi, axis=1), 1e-30)[:, None]
        else:
            M = self.Ms * self.m[c[:, 0], c[:, 1], c[:, 2]]
        return M * inside[:, None], inside


def make_mesh(spec, a):
    if spec.startswith("box:"):
        # the periodic box mesh of Firedrake (glueplex.PeriodicMesh does not accept the facet orientations of
        # fd.BoxMesh); is_periodic makes DemagField use the periodic solve (as for glueplex meshes)
        n = int(spec[4:])
        mesh = fd.PeriodicBoxMesh(n, n, n, a, a, a, reorder=False)
        mesh.is_periodic = True
        return mesh
    return PeriodicMesh(fd.Mesh(spec, reorder=False))


def cell_vertices(mesh):
    """Coordinates of the 4 vertices of each cell [ncell, 4, 3] (periodic meshes have a discontinuous coordinate
    field; the vertices of one cell are unwrapped)."""
    X = mesh.coordinates
    cmap = X.cell_node_map().values
    return X.dat.data_ro[cmap]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--body", required=True, help="vox, smooth or q<k> (the k-staircase, FD exact: X3)")
    ap.add_argument("--m", choices=["state", "interp", "100", "111"], default="state",
                    help="state: the coarse m (constant per FD cell); interp: a continuous m from it (trilinear in "
                         "the cell centres, unit length); 100, 111: uniform")
    ap.add_argument("--levels", type=int, default=2)
    ap.add_argument("--levels_cut", type=int, default=5, help="smooth body: sampling levels of the cells at the surface")
    ap.add_argument("--rtol", type=float, default=1e-10)
    ap.add_argument("--degree", type=int, default=1, help="degree of the potential (P1 or P2)")
    ap.add_argument("--mfamily", choices=["DG0", "DG1"], default=None,
                    help="M per cell: DG0 (mean; exact load for P1) or DG1 (L2 projection; exact load for P2); "
                         "default: DG<degree-1>")
    ap.add_argument("--out", default=None, help="append one JSON line with the result to this file")
    a = ap.parse_args()
    data = pathlib.Path(a.data)
    a.mfamily = a.mfamily or ("DG0" if a.degree == 1 else "DG1")
    t0 = time.time()
    src = Source(data, a.body, a.m)
    mesh = make_mesh(a.mesh, src.a)
    V0 = fd.FunctionSpace(mesh, "DG", 0)
    if a.mfamily == "DG0":
        W0 = fd.VectorFunctionSpace(mesh, "DG", 0)
    else:
        # nodes at the vertices (the default DG1 variant has inner nodes), checked against the coordinate field
        W0 = fd.FunctionSpace(mesh, fd.VectorElement(fd.FiniteElement("DG", mesh.ufl_cell(), 1, variant="equispaced")))
        xn = fd.Function(W0).interpolate(fd.SpatialCoordinate(mesh)).dat.data_ro[W0.cell_node_map().values]
        err = np.abs(xn - cell_vertices(mesh)).max()
        if err > 1e-9 * src.a:
            raise RuntimeError("DG1 nodes are not the cell vertices in their order (max distance %.3g m)" % err)
    verts = cell_vertices(mesh)
    ncell = len(verts)
    log("mesh %s: %d cells, %d vertices (%.0f s)" % (a.mesh, ncell, mesh.num_vertices(), time.time() - t0))

    t = time.time()
    nb = 1 if a.mfamily == "DG0" else 4
    coef = np.zeros((ncell, nb, 3))      # DG0: the mean of M; DG1: the L2 projection of M onto P1 in each cell
    frac = np.zeros(ncell)
    minv = np.linalg.inv((np.ones((4, 4)) + np.eye(4)) / 20.0)        # inverse P1 mass matrix of a unit-volume cell

    def sample(cells, lam):
        chunk = max(1, 4000000 // len(lam))
        for i in range(0, len(cells), chunk):
            c = cells[i:i + chunk]
            x = np.einsum("pk,ckd->cpd", lam, verts[c]).reshape(-1, 3)
            M, ins = src(x)
            M = M.reshape(len(c), len(lam), 3)
            if nb == 1:
                coef[c, 0] = M.mean(1)
            else:
                mom = np.einsum("pk,cpd->ckd", lam, M) / len(lam)             # mean of M λ_k
                coef[c] = np.einsum("jk,ckd->cjd", minv, mom)
            frac[c] = ins.reshape(len(c), len(lam)).mean(1)
    lam = ref_points(a.levels)
    sample(np.arange(ncell), lam)
    if a.body == "smooth" and a.levels_cut > a.levels:
        cut = np.flatnonzero(src.cut_cells(verts))                  # cells that the sphere surfaces may cut
        sample(cut, ref_points(a.levels_cut))
        log("%d cells near the sphere surfaces sampled with %d points" % (len(cut), 8 ** a.levels_cut))
    vol = np.abs(np.einsum("cd,cd->c", np.cross(verts[:, 1] - verts[:, 0], verts[:, 2] - verts[:, 0]),
                           verts[:, 3] - verts[:, 0])) / 6.0
    phi_mesh = float((frac * vol).sum() / vol.sum())
    log("sources: %d points per cell, magnetic volume fraction of the mesh %.6f, box volume error %.2e (%.0f s)"
        % (len(lam), phi_mesh, vol.sum() / src.a ** 3 - 1, time.time() - t))

    state = mp.State(mesh, default_degree=a.degree)
    m_fn = fd.Function(W0, name="m")
    dofs = W0.cell_node_map().values               # DG1: local nodes in the vertex order of the coordinate field
    for k in range(nb):
        m_fn.dat.data[dofs[:, k]] = coef[:, k] / src.Ms
    state.m = m_fn
    state.material.Ms = fd.Constant(src.Ms)
    PETSc.Options()["AirboxSolver_ksp_rtol"] = a.rtol
    demag = mp.DemagField()
    t = time.time()
    u = demag.u(state)          # potential (P1); DemagField.E needs m in the P1 space, so the energy is assembled here
    t_solve = time.time() - t
    E = fd.assemble(0.5 * MU0 * src.Ms * fd.inner(fd.grad(u), m_fn) * fd.dx)      # −½ μ0 ∫ M·H, H = −∇u
    E_grad = fd.assemble(0.5 * MU0 * fd.inner(fd.grad(u), fd.grad(u)) * fd.dx)    # ½ μ0 ∫ |∇u|²: equal if solved
    Kd = src.meta["Kd"]
    e = E / (Kd * src.a ** 3)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    res = {"data": str(data), "mesh": a.mesh, "body": a.body, "m": a.m, "levels": a.levels, "levels_cut": a.levels_cut if a.body == "smooth" else None, "degree": a.degree, "mfamily": a.mfamily, "cells": ncell,
           "vertices": mesh.num_vertices(), "E_KdVbox": e, "E_grad_rel": E_grad / E - 1, "phi_mesh": phi_mesh,
           "t_solve_s": t_solve,
           "t_total_s": time.time() - t0, "rss_GB": rss}
    if a.m not in ("state", "interp"):
        res["target"] = phi_mesh * (1 - phi_mesh) / 3
    log("E %.9f K_d V_box; solve %.0f s, total %.0f s, peak RSS %.2f GB" % (e, t_solve, time.time() - t0, rss))
    print(json.dumps(res))
    if a.out:
        with open(a.out, "a") as f:
            f.write(json.dumps(res) + "\n")


if __name__ == "__main__":
    main()
