#!/bin/bash
# box 11, PLAN 14b D7, cause of the G fail in q4_108 (amendment item 7; PROTOKOLL 3.41). The old-new difference of
# <Ms m_x> sits in the sample intervals 27 ... 29 and 99 ... 102. The chain c63e862 -> ABP has 6 steps; the runs below
# put the intermediate trees of the chain in place of ABP, with the same arguments and the same demag kernel (the
# kernel files of all trees have the same md5): B_d53P (PLAN 10, B5, B4 with the pair sum; no B1, no A), B_det (B_df0P:
# adds B1, the approved model change) and BA4P (adds A4). Each is judged by G against the old runs and the null runs of
# q4_108 as ABP. Information: two more runs of ABP in off108 (the run off108_new_r2 is 16 % slower than r1 and r3).
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"; eval "$(grep '^ARGS_A=' /root/s/k1.sh)"
A="${ARGS_A/--precision_cycles double /}"; A="${A/--init_from \/root\/states\/K1_2_sat\/init.pt/}"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
W=/root/d7c; R=/root/d7cres; D7=/root/d7; mkdir -p $W $R
run() {  # NAME TREE KC ARGS
  n=$1; t=$2; [ -f $W/$n.json ] && return
  [ -d /root/kc/c_$3 ] || cp -a /root/kc/$3 /root/kc/c_$3
  echo "$(date -u +%T) start D7c $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d7c_$n PYTHONPATH=/root/$t python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/c_$3 --json $W/$n.json -- $4 --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D7c $n exit $rc" >> $P
}
for t in B_d53P B_det BA4P; do run q4_108_$t $t ABP_108 "$B --init_from /root/states/Bprep/init.pt"; done
for r in 4 5; do run off108_new_r$r ABP ABP_108 "$A --init_from /root/states/K1_2_sat/init.pt"; done
for t in B_d53P B_det BA4P; do
  python /root/s/k_eval.py --g_old $D7/q4_108_old_r1.json $D7/q4_108_old_r2.json $D7/q4_108_old_r3.json \
    --g_new $W/q4_108_$t.json --g_null $D7/q4_108_null_reverse.json $D7/q4_108_null_rotate.json \
    --json $R/g_q4_108_$t.json > $R/g_q4_108_$t.txt 2>&1
  rc=$?; echo "$(date -u +%T) D7c G q4_108 $t exit $rc" >> $P
done
echo "$(date -u +%T) DONE d7c" >> $P
