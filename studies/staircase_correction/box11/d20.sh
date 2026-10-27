#!/bin/bash
# box 11, PLAN 16 (PROTOKOLL 3.46): the checks of AS5 (AS3 + G + Q) on the GPU, one job at a time.
# E1 (class E, bitwise): T1 and T2 of AS3 (again) and AS5 on the original B and A states, with copies of one inductor
#    cache; r16_eval.py --same on the T1 state files and on the T2 files (also against the T2 file of AS3 of PLAN 15);
#    speed runs of AS5 as W5 (r1) in off108, off156, q4_108, q4_156, q8_108: the trace against the 3 AS3 runs of W5.
# C1 (class R, the decision on new data): T1 and T2 of ABP and AS5 on the states of refB.pt and refA.pt with m + 1e-7 r
#    (perturb_ref.py, seed 3); k_eval.py --t1 (R1) and r16_eval.py --t2 (R2, R3).
# The T1 folders are deleted after their evaluation (about 7 GB each; as d12.sh). d21.sh (E2, CPU) runs at the same
# time; d22.sh (WQ) waits for both.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_A=' /root/s/k1.sh)"; eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
A32="${ARGS_A/--precision_cycles double /}"; B32="${ARGS_B/--precision_cycles double /}"
A16="${A32/--init_from \/root\/states\/K1_2_sat\/init.pt/}"; B16="${B32/--init_from \/root\/states\/Bprep\/init.pt/}"
args() { if [ $1 = A ]; then [ $2 = 64 ] && echo "$ARGS_A" || echo "$A32"; else [ $2 = 64 ] && echo "$ARGS_B" || echo "$B32"; fi; }
W=/root/p16; RES=$W/res; mkdir -p $W/ref $W/t1 $W/t2 $RES
log() { echo "$(date -u +%T) $*" >> $P; }
ic() {  # NAME GEO PREC: a copy of the base inductor cache of W2 of this geometry and precision
  d=/root/ic/p16_$1_$2_$3; [ -d $d ] || cp -a /root/ic/w2_$2/$3 $d; echo $d; }
t1() {  # TAG TREE GEO REF
  o=$W/t1/$1_$2_$3; [ -f $o/state_240.pt ] && return
  log "start P16 T1 $1 $2 $3"; cd /root/$2/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=$(ic $1_$2 $3 64) PYTHONPATH=/root/$2 python /root/s/s1_states.py --ref $4 --out $o -- $(args $3 64) --b_list 0.006 --out $W/run_t1_$1_$2_$3 > $o.log 2>&1
  rc=$?; log "end P16 T1 $1 $2 $3 exit $rc"; }
t2() {  # TAG TREE GEO REF
  o=$W/t2/$1_$2_$3.pt; [ -f $o ] && return
  log "start P16 T2 $1 $2 $3"; cd /root/$2/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=$(ic $1_$2 $3 32) PYTHONPATH=/root/$2 python /root/s/s1_states.py --ref $4 --out $o --against $W/t1/$1_$2_$3 --reps 3 -- $(args $3 32) --b_list 0.006 --out $W/run_t2_$1_$2_$3 > ${o%.pt}.log 2>&1
  rc=$?; log "end P16 T2 $1 $2 $3 exit $rc"; }
ev() {  # NAME COMMAND...: an evaluation with its exit code in the progress file
  n=$1; shift; "$@" > $RES/$n.txt 2>&1; rc=$?; log "P16 eval $n exit $rc"; }
# new states (seed 3)
for g in A B; do
  [ -f $W/ref/ref${g}_p3.pt ] || { python /root/s/perturb_ref.py /root/k1/$g/ref$g.pt $W/ref/ref${g}_p3.pt --seed 3 --eps 1e-7 > $W/ref/perturb_$g.txt 2>&1; log "P16 perturb_ref $g exit $?"; }
done
# E1: original states, AS3 again and AS5
for g in B A; do
  R=/root/k1/$g/ref$g.pt
  t1 e1 AS3 $g $R; t1 e1 AS5 $g $R
  ev e1_same_t1_$g python /root/s/r16_eval.py --same $W/t1/e1_AS3_$g $W/t1/e1_AS5_$g
  t2 e1 AS3 $g $R; t2 e1 AS5 $g $R
  ev e1_same_t2_$g python /root/s/r16_eval.py --same $W/t2/e1_AS3_$g.pt $W/t2/e1_AS5_$g.pt
  ev e1_same_t2_w2_$g python /root/s/r16_eval.py --same /root/w2/t2/$g/AS3.pt $W/t2/e1_AS5_$g.pt
  [ -f $RES/e1_same_t2_$g.txt ] && rm -rf $W/t1/e1_AS3_$g $W/t1/e1_AS5_$g
done
# E1: speed runs of AS5 (as d14.sh, r1), the trace against AS3 of W5
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
declare -A CONF NN ST
CONF[off108]="$A16"; NN[off108]=108; ST[off108]=K1_2_sat
CONF[off156]="$(N156 "$A16")"; NN[off156]=156; ST[off156]=K4_2_sat
CONF[q4_108]="$B16"; NN[q4_108]=108; ST[q4_108]=Bprep
CONF[q4_156]="$(N156 "$B16")"; NN[q4_156]=156; ST[q4_156]=T_q4_156
CONF[q8_108]="${B16/--staircase_q 4/--staircase_q 8}"; NN[q8_108]=108; ST[q8_108]=T_q8_108
for N in 108 156; do [ -d /root/kc/AS5_$N ] || cp -a /root/kc/AS3_$N /root/kc/AS5_$N; done
for c in q4_108 off108 q8_108 q4_156 off156; do
  n=${c}_AS5_r1; f=$W/$n.json
  if [ ! -f $f ]; then
    log "start P16 speed $n"; cd /root/AS5/studies/beta_bpeak
    TORCHINDUCTOR_CACHE_DIR=/root/ic/p16_$n PYTHONPATH=/root/AS5 python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/AS5_${NN[$c]} --json $f -- ${CONF[$c]} --init_from /root/states/${ST[$c]}/init.pt --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
    log "end P16 speed $n exit $?"
  fi
  ev e1_trace_$c python -c "
import json, sys
t = json.load(open('$f'))['trace']
o = [json.load(open('/root/w5/${c}_AS3_r%d.json' % r))['trace'] for r in (1, 2, 3)]
eq = [t == x for x in o]
print('trace $c AS5 r1 against AS3 r1, r2, r3 of W5 (126 rows: t, tries, <Ms m_x>, accepted):', eq)
sys.exit(0 if all(eq) and len(t) == 126 else 1)"
done
# C1: new states, ABP and AS5
for g in B A; do
  R=$W/ref/ref${g}_p3.pt
  t1 c1 ABP $g $R; t1 c1 AS5 $g $R
  ev c1_r1_$g python /root/s/k_eval.py --t1 R1_$g $W/t1/c1_ABP_$g $W/t1/c1_AS5_$g --json $RES/c1_r1_$g.json
  t2 c1 ABP $g $R; t2 c1 AS5 $g $R
  ev c1_r23_$g python /root/s/r16_eval.py --t2 AS5_$g $W/t2/c1_AS5_$g.pt $W/t2/c1_ABP_$g.pt --json $RES/c1_r23_$g.json
  [ -f $RES/c1_r23_$g.txt ] && rm -rf $W/t1/c1_ABP_$g $W/t1/c1_AS5_$g
done
# information: R2 and R3 on the PLAN 15 files (original B states, AS3 against ABP)
ev info_r23_w2_B python /root/s/r16_eval.py --t2 AS3_w2_B /root/w2/t2/B/AS3.pt /root/w2/t2/B/ABP.pt
log "DONE d20"
