#!/bin/bash
# box 11, PLAN 15 (PROTOKOLL 3.44), after d13.sh: (1) W5, the speed measurement as D7 (prof_step.py --skip 30 --plain
# 96 --window 0 --trace; one run on the GPU at a time; ABP and AS3 in turn; n 3 for each code), configurations off108,
# off156, q4_108, q4_156, q8_108, q8_156; the unperturbed runs are also realisation 1 of W4. (2) AS3 from the
# perturbed starts of d10a.sh (seeds 1 and 2) in off108, off156, q4_108 and q4_156 (realisations 2 and 3 of W4).
export CUDA_DEVICE=0
P=/root/progress.txt
while ! grep -q "DONE d13$" $P; do sleep 30; done
eval "$(grep '^ARGS_A=' /root/s/k1.sh)"; eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
A="${ARGS_A/--precision_cycles double /}"; A="${A/--init_from \/root\/states\/K1_2_sat\/init.pt/}"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
declare -A CONF NN ST
CONF[off108]="$A"; NN[off108]=108; ST[off108]=K1_2_sat
CONF[off156]="$(N156 "$A")"; NN[off156]=156; ST[off156]=K4_2_sat
CONF[q4_108]="$B"; NN[q4_108]=108; ST[q4_108]=Bprep
CONF[q4_156]="$(N156 "$B")"; NN[q4_156]=156; ST[q4_156]=T_q4_156
CONF[q8_108]="${B/--staircase_q 4/--staircase_q 8}"; NN[q8_108]=108; ST[q8_108]=T_q8_108
CONF[q8_156]="$(N156 "${B/--staircase_q 4/--staircase_q 8}")"; NN[q8_156]=156; ST[q8_156]=T_q8_156
for N in 108 156; do [ -d /root/kc/AS3_$N ] || cp -a /root/kc/ABP_$N /root/kc/AS3_$N; done
W=/root/w5; mkdir -p $W
run() {  # NAME TREE CONF START
  n=$1; t=$2; c=$3; [ -f $W/$n.json ] && return
  echo "$(date -u +%T) start D14 $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d14_$n PYTHONPATH=/root/$t python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/${t}_${NN[$c]} --json $W/$n.json -- ${CONF[$c]} --init_from $4 --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D14 $n exit $rc" >> $P
}
for r in 1 2 3; do
  for c in off108 off156 q4_108 q4_156 q8_108 q8_156; do
    for t in ABP AS3; do run ${c}_${t}_r$r $t $c /root/states/${ST[$c]}/init.pt; done
  done
done
for c in off108 off156 q4_108 q4_156; do for s in 1 2; do
  run ${c}_AS3_p$s AS3 $c /root/states/P${s}_${ST[$c]}/init.pt
done; done
echo "$(date -u +%T) DONE d14" >> $P
