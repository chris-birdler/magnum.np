#!/bin/bash
# box 11, PLAN 15 WFQ (PROTOKOLL 3.44): the set-up of ABP and of ASQ (FQ) bitwise compared with fq_check.py, at q 8
# for N 108 and N 156 (the configurations of the q 8 arms of D7) and at q 16 on the T3 geometry (d/l_ex 92, phi 0.55,
# N 48), on the GPU. Waits for the end of d10b.sh.
export CUDA_DEVICE=0
P=/root/progress.txt
while ! grep -q "DONE d10b$" $P; do sleep 30; done
W=/root/d10; mkdir -p $W
chk() {  # NAME CONFIG [OPTIONS]
  n=$1; [ -f $W/$n.json ] && return
  echo "$(date -u +%T) start D10c $n" >> $P
  cd /root/s
  PYTHONPATH=/root/ASQ:/root/ASQ/studies/beta_bpeak python /root/s/fq_check.py --old /root/ABP/magnumnp/field_terms/demag_staircase.py --new /root/ASQ/magnumnp/field_terms/demag_staircase.py --config $2 $3 --json $W/$n.json > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D10c $n exit $rc" >> $P
}
chk wfq_q8_108 /root/d7/run_q8_108_new_r1/config.json
chk wfq_q8_156 /root/d7/run_q8_156_new_r1/config.json
chk wfq_q16_T3 /root/d7/run_q4_108_new_r1/config.json "--q 16 --d_lex 92 --phi 0.55"
echo "$(date -u +%T) DONE d10c" >> $P
