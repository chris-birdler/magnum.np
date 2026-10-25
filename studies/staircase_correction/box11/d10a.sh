#!/bin/bash
# box 11, PLAN 15 step 1 (PROTOKOLL 3.44): runs that do not depend on AS. (1) The perturbed start files of W4
# (perturb_init.py, seeds 1 and 2, eps 1e-7) for the 4 configurations. (2) ABP from the perturbed starts (single
# precision, as the ABP arms of D7). (3) Dref of off108 and off156: 9e0ea7c in double precision; ABP in double
# precision in off108 and off156 (check of Dref). Each double run has its own kernel folder (as d7c3.sh).
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_A=' /root/s/k1.sh)"; eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
A2="${ARGS_A/--init_from \/root\/states\/K1_2_sat\/init.pt/}"; A="${A2/--precision_cycles double /}"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
declare -A CONF NN ST
CONF[off108]="$A"; NN[off108]=108; ST[off108]=K1_2_sat
CONF[off156]="$(N156 "$A")"; NN[off156]=156; ST[off156]=K4_2_sat
CONF[q4_108]="$B"; NN[q4_108]=108; ST[q4_108]=Bprep
CONF[q4_156]="$(N156 "$B")"; NN[q4_156]=156; ST[q4_156]=T_q4_156
W=/root/d10; mkdir -p $W
for c in off108 off156 q4_108 q4_156; do for s in 1 2; do
  d=/root/states/P${s}_${ST[$c]}; [ -f $d/init.pt ] && continue
  python /root/s/perturb_init.py /root/states/${ST[$c]}/init.pt $d/init.pt --seed $s --eps 1e-7 >> $W/perturb.txt 2>&1
  rc=$?; echo "$(date -u +%T) D10a perturb ${ST[$c]} seed $s exit $rc" >> $P
done; done
run() {  # NAME TREE KC ARGS
  n=$1; t=$2; k=$3; [ -f $W/$n.json ] && return
  echo "$(date -u +%T) start D10a $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d10_$n PYTHONPATH=/root/$t python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache $k --json $W/$n.json -- $4 --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D10a $n exit $rc" >> $P
}
for N in 108 156; do [ -d /root/kc/c10_ABP_$N ] || cp -a /root/kc/ABP_$N /root/kc/c10_ABP_$N; done
for c in off108 off156 q4_108 q4_156; do for s in 1 2; do
  run ${c}_ABP_p$s ABP /root/kc/c10_ABP_${NN[$c]} "${CONF[$c]} --init_from /root/states/P${s}_${ST[$c]}/init.pt"
done; done
OA108="$A2 --init_from /root/states/K1_2_sat/init.pt"; OA156="$(N156 "$A2") --init_from /root/states/K4_2_sat/init.pt"
for n in off108_M9e0_double off108_ABP_double off156_M9e0_double off156_ABP_double; do mkdir -p /root/kc/c_d_$n; done
run off108_M9e0_double M_9e0 /root/kc/c_d_off108_M9e0_double "$OA108"
run off108_ABP_double ABP /root/kc/c_d_off108_ABP_double "$OA108"
run off156_M9e0_double M_9e0 /root/kc/c_d_off156_M9e0_double "$OA156"
run off156_ABP_double ABP /root/kc/c_d_off156_ABP_double "$OA156"
echo "$(date -u +%T) DONE d10a" >> $P
