#!/bin/bash
# PLAN 4.21 K6, voxel regression of the evaluation patch: the post-processing scripts of an old commit (a temporary
# worktree) and of this worktree on two copies of the same finished voxel run (with accumulators and snapshots); all
# outputs must be equal (text of the JSON and md files, arrays of the npz files). Usage:
#   check421_voxel_regression.sh <old commit> <run folder> <python> [<scratch folder>]
# dissmap_batch.py runs also on a copy with only the 4 phases 90, 180, 270, 360 (dissmap4.json).
# Result 2026-11-06: e823a8e against the patch, a K5 voxel run (d/l_ex 48, 2 amplitudes): REGRESSION PASS (all
# outputs bitwise equal; dissmap4.json: R within 1e-7, the float32 rounding of the mean, see below).
OLDC=$1; RUN=$2; V=$3; S=${4:-$(mktemp -d)}
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
NEW=$REPO/studies/beta_bpeak
git -C $REPO worktree add -q --detach $S/wt_old $OLDC || exit 1
OLD=$S/wt_old/studies/beta_bpeak
rm -rf $S/reg; mkdir -p $S/reg/old $S/reg/new
cp -a $RUN $S/reg/old/K5v; cp -a $RUN $S/reg/new/K5v
for side in old new; do
  [ $side = old ] && C=$OLD || C=$NEW
  [ $side = old ] && PP=$S/wt_old || PP=$REPO
  d=$S/reg/$side/K5v
  cd $C
  for sc in acc_reduce.py phasemap_batch.py phase_extra.py domains_batch.py dissmap_batch.py; do
    OMP_NUM_THREADS=4 CUDA_DEVICE=-1 PYTHONPATH=$PP $V $sc $d > $S/reg/$side.$sc.log 2>&1; echo "$side $sc exit $?"
  done
  OMP_NUM_THREADS=4 CUDA_DEVICE=-1 PYTHONPATH=$PP $V eval_hotspots.py --json $d/hotspots.json $d > $d/hotspots.md 2> $S/reg/$side.hs.err; echo "$side hotspots exit $?"
  OMP_NUM_THREADS=4 CUDA_DEVICE=-1 PYTHONPATH=$PP $V eval_hotcores.py --json $d/hotcores.json $d > $d/hotcores.md 2> $S/reg/$side.hc.err; echo "$side hotcores exit $?"
  OMP_NUM_THREADS=4 CUDA_DEVICE=-1 PYTHONPATH=$PP $V census_bloch.py --out $d/census_bloch.json $d > $S/reg/$side.cb.log 2>&1; echo "$side census exit $?"
  # dissmap_batch.py needs exactly the 4 phases 90, 180, 270, 360 of a cycle: a copy with only these snapshots
  d4=$S/reg/$side/K5v4; rm -rf $d4; cp -a $RUN $d4
  find $d4 -name 'm_*_ph*.vti' ! -name '*_ph090.vti' ! -name '*_ph180.vti' ! -name '*_ph270.vti' ! -name '*_ph360.vti' -delete
  OMP_NUM_THREADS=4 CUDA_DEVICE=-1 PYTHONPATH=$PP $V dissmap_batch.py $d4 > $S/reg/$side.dm4.log 2>&1; echo "$side dissmap (4 phases) exit $?"
  cp $d4/dissmap.json $d/dissmap4.json
done
$V - $S/reg <<'PY'
import json, pathlib, sys, numpy as np
S = pathlib.Path(sys.argv[1])
o, n = S / "old" / "K5v", S / "new" / "K5v"
bad = 0
for f in sorted(o.glob("*.json")) + sorted(o.glob("*.md")):
    if f.name in ("config.json", "summary.json"):
        continue
    a = f.read_text().replace(str(o), "RUN").replace("reg/old", "reg/X")
    b = (n / f.name).read_text().replace(str(n), "RUN").replace("reg/new", "reg/X") if (n / f.name).exists() else None
    same = a == b
    note = ""
    if not same and f.name.startswith("dissmap") and b is not None:
        # dissmap_batch.py: R = (float32 mean of p_f) / p_dis. The float32 sum of numpy depends on the memory
        # alignment of the array (2026-11-06: two equal arrays at offsets 16 and 32 of 64 gave means 1 rounding
        # step apart), so R may differ by about 1e-7 also for equal inputs. Numbers within 1e-6 relative pass.
        def flat(x, out):
            if isinstance(x, dict):
                for k in sorted(x):
                    flat(x[k], out)
            elif isinstance(x, list):
                for v in x:
                    flat(v, out)
            else:
                out.append(x)
            return out
        A, Bv = flat(json.loads(a), []), flat(json.loads(b), [])
        same = len(A) == len(Bv) and all(x == y or (isinstance(x, float) and isinstance(y, float) and
                                                      abs(x - y) <= 1e-6 * max(abs(x), abs(y))) for x, y in zip(A, Bv))
        note = " (within the float32 rounding of the mean, 1e-6 relative)" if same else ""
    print("%-22s %s%s" % (f.name, "equal" if same else "DIFFERENT", note))
    bad += not same
for f in sorted(o.glob("*.npz")):
    if f.name.startswith("acc_") and f.name != "acc_slices.npz" or f.name == "geometry.npz":
        continue
    za, zb = np.load(f), np.load(n / f.name)
    keys = sorted(set(za.files) | set(zb.files))
    diff = [k for k in keys if k not in za.files or k not in zb.files or not np.array_equal(za[k], zb[k], equal_nan=True)]
    print("%-22s %s" % (f.name, "equal (%d arrays)" % len(keys) if not diff else "DIFFERENT %s" % diff[:5]))
    bad += bool(diff)
print("REGRESSION", "PASS" if bad == 0 else "FAIL (%d)" % bad)
sys.exit(1 if bad else 0)
PY
rc=$?
git -C $REPO worktree remove --force $S/wt_old
exit $rc
