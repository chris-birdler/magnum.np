#!/bin/bash
# box 11, PLAN 14b D6: T1, T2 and T6 of the B chain with the P trees (states of refB.pt, criteria of PLAN 14a). The T2
# files of c63e862 (B_c63, B_c63n in /root/t2/B, from t12.sh) stay valid; T1 of c63e862 runs again (its folder was
# deleted after use). At most 3 double-precision folders on the disk at a time.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
B32="${ARGS_B/--precision_cycles double /}"
T1=/root/t1/BP; T2=/root/t2/B; RES=/root/tres_p
mkdir -p $T1 $RES
r64() {  # NAME TREE
  n=$1; t=$2; [ -f $T1/$n/state_240.pt ] && return
  echo "$(date -u +%T) start T1 BP $n" >> $P; cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/t1_BP_$n PYTHONPATH=/root/$t python /root/s/s1_states.py --ref /root/k1/B/refB.pt --out $T1/$n -- $ARGS_B --b_list 0.006 --out /root/t1/run_BP_$n > $T1/$n.log 2>&1
  echo "$(date -u +%T) end T1 BP $n exit $?" >> $P; }
r32() {  # NAME TREE
  n=$1; t=$2; [ -f $T2/$n.pt ] && return
  echo "$(date -u +%T) start T2 BP $n" >> $P; cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/t2_BP_$n PYTHONPATH=/root/$t python /root/s/s1_states.py --ref /root/k1/B/refB.pt --out $T2/$n.pt --against $T1/$n --reps 3 -- $B32 --b_list 0.006 --out /root/t2/run_BP_$n > $T2/$n.log 2>&1
  echo "$(date -u +%T) end T2 BP $n exit $?" >> $P; }
cmp1() {  # NAME OLD NEW
  python /root/s/k_eval.py --t1 $1 $T1/$2 $T1/$3 --json $RES/t1_$1.json > $RES/t1_$1.txt 2>&1
  echo "$(date -u +%T) T1 P $1 exit $?" >> $P; }
r64 B_c63 B_c63
r64 B_01eP B_01eP; cmp1 PLAN10 B_c63 B_01eP; r32 B_01eP B_01eP; rm -rf $T1/B_c63
r64 B_c56P B_c56P; cmp1 B5 B_01eP B_c56P; r32 B_c56P B_c56P; rm -rf $T1/B_01eP
r64 B_d53P B_d53P; cmp1 B4 B_c56P B_d53P; r32 B_d53P B_d53P; rm -rf $T1/B_c56P
r64 B_df0P B_det; python /root/s/k_eval.py --t6 $T1/B_d53P $T1/B_df0P --json $RES/t6.json > $RES/t6.txt 2>&1
echo "$(date -u +%T) T6 P exit $?" >> $P; r32 B_df0P B_det; rm -rf $T1/B_d53P
r64 BA4P BA4P; cmp1 A4 B_df0P BA4P; r32 BA4P BA4P; rm -rf $T1/B_df0P
r64 ABP_on ABP; cmp1 AB_on BA4P ABP_on; r32 ABP_on ABP; rm -rf $T1/BA4P $T1/ABP_on
# T2 evaluation: NAME NEW OLD (chain start c63e862 and its null run)
t2e() { python /root/s/k_eval.py --t2 $1 $T2/$2.pt $T2/$3.pt $T2/B_c63.pt $T2/B_c63n.pt --json $RES/t2_$1.json > $RES/t2_$1.txt 2>&1; echo "$(date -u +%T) T2 P $1 exit $?" >> $P; }
t2e PLAN10 B_01eP B_c63; t2e B5 B_c56P B_01eP; t2e B4 B_d53P B_c56P; t2e B1 B_df0P B_d53P; t2e A4 BA4P B_df0P
t2e AB_on ABP_on BA4P
touch $RES/T12P_DONE
echo "$(date -u +%T) DONE t12p" >> $P
