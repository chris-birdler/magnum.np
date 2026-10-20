#!/bin/bash
# FE comparison, preparation (PROTOKOLL 3.24): night run on the local machine (Chris, go 2026-10-20).
# Host side: the venv of magnum.np (VENV, with torch and scipy); FE side: docker image magnum.pi.dev-3.13.3 with the
# private magnum.pi clone mounted (never published). Results: fe_data/fe_night.jsonl, results/fe_static.{md,json}.
set -u
VENV=${VENV:?path of the python venv}
CACHE=${CACHE:?folder of the FD kernel cache}
MPI=${MPI:-$HOME/Documents/code/magnum.pi}
W=$(cd "$(dirname "$0")/.." && pwd)
ROOT=$(cd "$W/../.." && pwd)
OUT=fe_data/fe_night.jsonl
STATES="T3Pc_s941:fe_data/src_states/T3Pc_s941 T3Pc_s942:fe_data/src_states/T3Pc_s942 T3Pc_s943:fe_data/src_states/T3Pc_s943 T3Pc_s944:fe_data/src_states/T3Pc_s944 X3_s951:remote_20261020_box3/runs/X3_s951 X3_s952:remote_20261020_box3/runs/X3_s952"
fe() { docker run --rm --memory 12g -v "$MPI":/home/firedrake/magnum.pi -v "$W":/work magnum.pi.dev-3.13.3 \
       bash -c "cd /work && python3 fe/fe_demag.py --out $OUT $*" 2>&1 | grep -E "^[0-9:]{8} |^\{|Error|rror|Killed"; }
cd "$W"
echo "== FD exports $(date -u +%H:%M)"
for s in $STATES; do
  PYTHONPATH="$ROOT" "$VENV/bin/python" fe/fe_states.py --run "${s#*:}" --out "fe_data/${s%%:*}" --q 2 4 --cache "$CACHE" \
    2>&1 | grep -E "^[0-9:]{8} |Error|rror"
done
echo "== FE 48^3, P2 $(date -u +%H:%M)"
for s in $STATES; do
  n=${s%%:*}
  fe --data fe_data/$n --mesh box:48 --body vox --m state --levels 0 --degree 2
  fe --data fe_data/$n --mesh box:48 --body q2 --m state --levels 3 --degree 2
  fe --data fe_data/$n --mesh box:48 --body smooth --m state --levels 0 --levels_cut 5 --degree 2
done
echo "== FE 48^3, continuous m $(date -u +%H:%M)"
for s in $STATES; do
  fe --data fe_data/${s%%:*} --mesh box:48 --body smooth --m interp --levels 2 --levels_cut 4 --degree 2
done
echo "== FE 64^3 (not conforming to the FD cells) $(date -u +%H:%M)"
for n in X3_s951 T3Pc_s941; do
  fe --data fe_data/$n --mesh box:64 --body vox --m state --levels 3 --degree 2
  fe --data fe_data/$n --mesh box:64 --body q2 --m state --levels 3 --degree 2
  fe --data fe_data/$n --mesh box:64 --body smooth --m state --levels 3 --levels_cut 5 --degree 2
done
echo "== evaluation $(date -u +%H:%M)"
"$VENV/bin/python" fe/fe_eval.py --results $OUT --out results/fe_static
echo "== NIGHT DONE $(date -u +%H:%M)"
