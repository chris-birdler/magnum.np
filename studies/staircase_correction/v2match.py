# Check of the V2 owner map for one voxel sphere (cell-centre rule) at q = 2:
# can the orphan fine cells of the volume-matched fine sphere be matched to the lost fine cells
# of the voxel cubes within a transport radius? Where do the unmatched orphans sit?
import sys, numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
R = float(sys.argv[1]); q = 2; radii = [float(v) for v in sys.argv[2].split(",")]
n = int(np.ceil(R)) + 4
c = np.arange(-n, n) + 0.5
X, Y, Z = np.meshgrid(c, c, c, indexing="ij"); vox = X**2 + Y**2 + Z**2 <= R**2; Nv = int(vox.sum())
f = (np.arange(-n*q, n*q) + 0.5) / q
FX, FY, FZ = np.meshgrid(f, f, f, indexing="ij"); r2 = FX**2 + FY**2 + FZ**2
vals = np.unique(r2[r2 <= (R + 2)**2]); cnt = np.searchsorted(np.sort(r2.ravel()), vals, side="right")
k = np.argmin(np.abs(cnt - q**3 * Nv)); body = r2 <= vals[k]
print("R %.3f cells: voxel cells %d, fine body %d, q^3 Nv %d, delta %d" % (R, Nv, body.sum(), q**3*Nv, body.sum()-q**3*Nv))
parent = np.repeat(np.repeat(np.repeat(vox, q, 0), q, 1), q, 2)
orph = body & ~parent; lost = parent & ~body
lost_c = lost.reshape(2*n, q, 2*n, q, 2*n, q).sum(axis=(1, 3, 5))
recv = np.argwhere(lost_c > 0); cap = lost_c[lost_c > 0].astype(int)
opos = (np.argwhere(orph) + 0.5) / q - n; rpos = recv + 0.5 - n
print("orphans %d, lost slots %d, receiving cells %d" % (len(opos), cap.sum(), len(rpos)), flush=True)
first = np.concatenate([[0], np.cumsum(cap)[:-1]]); tr = cKDTree(rpos)
u = opos / np.linalg.norm(opos, axis=1)[:, None]; pol = np.abs(u).max(axis=1)   # 1 at the 6 poles, 0.577 at (111)
for rad in radii:
    rows, cols = [], []
    for o, lst in enumerate(tr.query_ball_point(opos, rad)):
        for j in lst:
            rows.extend([o] * cap[j]); cols.extend(range(first[j], first[j] + cap[j]))
    A = csr_matrix((np.ones(len(rows), np.int8), (rows, cols)), shape=(len(opos), int(cap.sum())))
    mt = maximum_bipartite_matching(A, perm_type="column")
    un = mt < 0
    h = np.histogram(pol[un], bins=[0.57, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0001])[0] if un.any() else []
    print("radius %.1f dx: unmatched %d; max|n_i| of unmatched in bins [.57,.7,.8,.9,.95,.99,1]: %s" % (rad, int(un.sum()), list(h)), flush=True)
