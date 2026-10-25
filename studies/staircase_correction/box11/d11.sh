#!/bin/bash
# box 11, PLAN 15 W2 (PROTOKOLL 3.44): T1 and T2 (the criteria of 14a) of the chain ABP -> ASQ -> AS1 -> AS2 -> AS3 on
# the 16 B states of refB.pt; for F2 and F3 also on the A states (refA.pt). The two trees of a pair run with copies of
# one inductor cache (X1: then block_sum gets the same kernel configuration); the chosen configurations of the stage-1
# kernel of block_sum are recorded. Usage: d11.sh PAIRNAME OLDTREE NEWTREE GEO(A|B)
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_A=' /root/s/k1.sh)"; eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
A32="${ARGS_A/--precision_cycles double /}"; B32="${ARGS_B/--precision_cycles double /}"
pair=$1; old=$2; new=$3; g=$4
ref() { [ $1 = A ] && echo /root/k1/A/refA.pt || echo /root/k1/B/refB.pt; }
args() { if [ $1 = A ]; then [ $2 = 64 ] && echo "$ARGS_A" || echo "$A32"; else [ $2 = 64 ] && echo "$ARGS_B" || echo "$B32"; fi; }
T1=/root/w2/t1/$g; T2=/root/w2/t2/$g; RES=/root/w2/res; IC=/root/ic/w2_$g; mkdir -p $T1 $T2 $RES
r64() {  # TREE
  t=$1; [ -f $T1/$t/state_240.pt ] && return
  [ -d $IC/64_$t ] || { [ -d $IC/64 ] && cp -a $IC/64 $IC/64_$t || mkdir -p $IC/64_$t; }
  echo "$(date -u +%T) start W2 T1 $g $t" >> $P; cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=$IC/64_$t PYTHONPATH=/root/$t python /root/s/s1_states.py --ref $(ref $g) --out $T1/$t -- $(args $g 64) --b_list 0.006 --out /root/w2/run_${g}_64_$t > $T1/$t.log 2>&1
  rc=$?; echo "$(date -u +%T) end W2 T1 $g $t exit $rc" >> $P
  [ -d $IC/64 ] || cp -a $IC/64_$t $IC/64; }
r32() {  # TREE
  t=$1; [ -f $T2/$t.pt ] && return
  [ -d $IC/32_$t ] || { [ -d $IC/32 ] && cp -a $IC/32 $IC/32_$t || mkdir -p $IC/32_$t; }
  echo "$(date -u +%T) start W2 T2 $g $t" >> $P; cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=$IC/32_$t PYTHONPATH=/root/$t python /root/s/s1_states.py --ref $(ref $g) --out $T2/$t.pt --against $T1/$t --reps 3 -- $(args $g 32) --b_list 0.006 --out /root/w2/run_${g}_32_$t > $T2/$t.log 2>&1
  rc=$?; echo "$(date -u +%T) end W2 T2 $g $t exit $rc" >> $P
  [ -d $IC/32 ] || cp -a $IC/32_$t $IC/32; }
cfg() {  # the configurations of the stage-1 kernel of block_sum (triton_poi_fused_stack_0) in the caches of a tree
  for d in $IC/64_$1 $IC/32_$1; do
    for f in $(grep -rl "def triton_poi_fused_stack_0" $d --include=*.py 2>/dev/null); do
      echo "$1 $(basename $d): $(cat $(dirname $f)/*.best_config 2>/dev/null | tr '\n' ' ')"; done; done; }
r64 $old; r64 $new
python /root/s/k_eval.py --t1 $pair $T1/$old $T1/$new --json $RES/t1_$pair.json > $RES/t1_$pair.txt 2>&1
rc=$?; echo "$(date -u +%T) W2 T1 $pair exit $rc" >> $P
r32 $old; r32 $new
# chain start and its null run as in 14a/14b: c63e862 (B) and 9b15ea9 (A), files of t12.sh
if [ $g = B ]; then ST=/root/t2/B/B_c63.pt; NL=/root/t2/B/B_c63n.pt; else ST=/root/t2/A/A_9b1.pt; NL=/root/t2/A/A_9b1n.pt; fi
python /root/s/k_eval.py --t2 $pair $T2/$new.pt $T2/$old.pt $ST $NL --json $RES/t2_$pair.json > $RES/t2_$pair.txt 2>&1
rc=$?; echo "$(date -u +%T) W2 T2 $pair exit $rc" >> $P
{ cfg $old; cfg $new; } > $RES/cfg_$pair.txt
echo "$(date -u +%T) DONE d11 $pair" >> $P
