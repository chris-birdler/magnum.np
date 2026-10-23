#!/bin/bash
# box 11, PLAN 14a T1 (operators in double precision), T2 (single-precision kernels), T6 (B1); chains with at most 3
# double-precision folders on the disk at a time. Needs the K1 references (/root/k1/A/refA.pt, /root/k1/B/refB.pt).
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_A=' /root/s/k1.sh)"; eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
A32="${ARGS_A/--precision_cycles double /}"; B32="${ARGS_B/--precision_cycles double /}"
mkdir -p /root/t1/A /root/t1/B /root/t2/A /root/t2/B /root/tres
ref() { [ $1 = A ] && echo /root/k1/A/refA.pt || echo /root/k1/B/refB.pt; }
args() { if [ $1 = A ]; then [ $2 = 64 ] && echo "$ARGS_A" || echo "$A32"; else [ $2 = 64 ] && echo "$ARGS_B" || echo "$B32"; fi; }
r64() {  # NAME TREE GEO [--null]
  n=$1; t=$2; g=$3; shift 3; [ -f /root/t1/$g/$n/state_240.pt ] && return
  echo "$(date -u +%T) start T1 $g $n" >> $P; cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/t1_${g}_$n PYTHONPATH=/root/$t python /root/s/s1_states.py --ref $(ref $g) --out /root/t1/$g/$n "$@" -- $(args $g 64) --b_list 0.006 --out /root/t1/run_${g}_$n > /root/t1/$g/$n.log 2>&1
  echo "$(date -u +%T) end T1 $g $n exit $?" >> $P; }
r32() {  # NAME TREE GEO AGAINST [--null]
  n=$1; t=$2; g=$3; ag=$4; shift 4; [ -f /root/t2/$g/$n.pt ] && return
  echo "$(date -u +%T) start T2 $g $n" >> $P; cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/t2_${g}_$n PYTHONPATH=/root/$t python /root/s/s1_states.py --ref $(ref $g) --out /root/t2/$g/$n.pt --against /root/t1/$g/$ag --reps 3 "$@" -- $(args $g 32) --b_list 0.006 --out /root/t2/run_${g}_$n > /root/t2/$g/$n.log 2>&1
  echo "$(date -u +%T) end T2 $g $n exit $?" >> $P; }
cmp1() {  # NAME GEO OLD NEW
  python /root/s/k_eval.py --t1 $1 /root/t1/$2/$3 /root/t1/$2/$4 --json /root/tres/t1_$1.json > /root/tres/t1_$1.txt 2>&1
  echo "$(date -u +%T) T1 $1 exit $?" >> $P; }
# A chain: 9b15ea9 -> 31aa455 -> A1 fixed -> 3d9de71 -> AB (correction off)
r64 A_9b1 A_9b1 A; r64 A_9b1n A_9b1 A --null; cmp1 A_null A A_9b1 A_9b1n
r32 A_9b1 A_9b1 A A_9b1; r32 A_9b1n A_9b1 A A_9b1 --null; rm -rf /root/t1/A/A_9b1n
r64 A_31a A_31a A; cmp1 A2 A A_9b1 A_31a; r32 A_31a A_31a A A_31a
r64 A_a1f A_a1f A; cmp1 A1 A A_31a A_a1f; r32 A_a1f A_a1f A A_a1f; rm -rf /root/t1/A/A_31a
r64 A_3d9 A_3d9 A; cmp1 A3 A A_a1f A_3d9; r32 A_3d9 A_3d9 A A_3d9; rm -rf /root/t1/A/A_a1f
r64 AB_off AB A; cmp1 AB_off A A_3d9 AB_off; r32 AB_off AB A AB_off; rm -rf /root/t1/A/A_3d9 /root/t1/A/AB_off /root/t1/A/A_9b1
# B chain: c63e862 -> 01ed90c -> c56c0a5 -> d538bde -> df0f811 (B1: T6) -> BA4 (A4) -> AB (correction on)
r64 B_c63 B_c63 B; r64 B_c63n B_c63 B --null; cmp1 B_null B B_c63 B_c63n
r32 B_c63 B_c63 B B_c63; r32 B_c63n B_c63 B B_c63 --null; rm -rf /root/t1/B/B_c63n
r64 B_01e B_01e B; cmp1 PLAN10 B B_c63 B_01e; r32 B_01e B_01e B B_01e; rm -rf /root/t1/B/B_c63
r64 B_c56 B_c56 B; cmp1 B5 B B_01e B_c56; r32 B_c56 B_c56 B B_c56; rm -rf /root/t1/B/B_01e
r64 B_d53 B_d53 B; cmp1 B4 B B_c56 B_d53; r32 B_d53 B_d53 B B_d53; rm -rf /root/t1/B/B_c56
r64 B_df0 B_df0 B; python /root/s/k_eval.py --t6 /root/t1/B/B_d53 /root/t1/B/B_df0 --json /root/tres/t6.json > /root/tres/t6.txt 2>&1
echo "$(date -u +%T) T6 exit $?" >> $P; r32 B_df0 B_df0 B B_df0; rm -rf /root/t1/B/B_d53
r64 BA4 BA4 B; cmp1 A4 B B_df0 BA4; r32 BA4 BA4 B BA4; rm -rf /root/t1/B/B_df0
r64 AB_on AB B; cmp1 AB_on B BA4 AB_on; r32 AB_on AB B AB_on; rm -rf /root/t1/B/BA4 /root/t1/B/AB_on
# T2 evaluation: NAME NEW OLD START NULL
t2e() { python /root/s/k_eval.py --t2 $1 /root/t2/$2/$3.pt /root/t2/$2/$4.pt /root/t2/$2/$5.pt /root/t2/$2/$6.pt --json /root/tres/t2_$1.json > /root/tres/t2_$1.txt 2>&1; echo "$(date -u +%T) T2 $1 exit $?" >> $P; }
t2e A2 A A_31a A_9b1 A_9b1 A_9b1n; t2e A1 A A_a1f A_31a A_9b1 A_9b1n; t2e A3 A A_3d9 A_a1f A_9b1 A_9b1n
t2e AB_off A AB_off A_3d9 A_9b1 A_9b1n
t2e PLAN10 B B_01e B_c63 B_c63 B_c63n; t2e B5 B B_c56 B_01e B_c63 B_c63n; t2e B4 B B_d53 B_c56 B_c63 B_c63n
t2e B1 B B_df0 B_d53 B_c63 B_c63n; t2e A4 B BA4 B_df0 B_c63 B_c63n; t2e AB_on B AB_on BA4 B_c63 B_c63n
touch /root/tres/T12_DONE
echo "$(date -u +%T) DONE t12" >> $P
