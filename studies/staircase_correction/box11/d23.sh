#!/bin/bash
# box 11, PLAN 16 (PROTOKOLL 3.46). (1) At once (CPU, small): calibration of R2 and R3 (r16_eval.py --t2) on the T2
# files of the A states of the PLAN 14 chain (start 9b15ea9, its null run, A2, A1, A3, AB off) and of PLAN 15 (AS1 ...
# AS3). (2) After d20.sh, on the GPU (E2 and L2; replaces d21.sh, whose CPU run was too slow under the CPU quota of the
# container): fq_check.py AS3 against AS5 (md5 of the centroid offsets, B, V, S, T, rows, _d, _v and the small-budget
# offsets), T3 geometry (d/l_ex 96, phi 0.55, N 48) at q 2, 4, 8, 16; N 108 at q 4, 8, 16; N 156 at q 4, 8.
export CUDA_DEVICE=0
P=/root/progress.txt
W=/root/p16; mkdir -p $W/cal $W/e2
T=/root/t2/A; V=/root/w2/t2/A
cal() {  # OLDFILE NEWFILE NAME
  python /root/s/r16_eval.py --t2 $3 $2 $1 > $W/cal/$3.txt 2>&1; echo "$(date -u +%T) P16 cal $3 exit $?" >> $P; }
cal $T/A_9b1.pt $T/A_9b1n.pt A_9b1-9b1n
cal $T/A_9b1n.pt $T/A_9b1.pt A_9b1n-9b1
for t in A_31a A_a1f A_3d9 AB_off; do cal $T/A_9b1.pt $T/$t.pt A_9b1-$t; done
cal $T/AB_off.pt $V/AS1.pt A_ABoff-AS1
cal $V/AS1.pt $V/AS2.pt A_AS1-AS2
cal $V/AS2.pt $V/AS3.pt A_AS2-AS3
cal $T/AB_off.pt $V/AS3.pt A_ABoff-AS3
cal $T/AB_off.pt $T/A_9b1.pt A_ABoff-9b1
while ! grep -q "DONE d20$" $P; do sleep 30; done
OLD=/root/AS3/magnumnp/field_terms/demag_staircase.py; NEW=/root/AS5/magnumnp/field_terms/demag_staircase.py
C108=/root/w5/run_q8_108_AS3_r1/config.json; C156=/root/w5/run_q8_156_AS3_r1/config.json
fq() {  # NAME CONFIG ARGS...
  n=$1; c=$2; shift 2; [ -f $W/e2/$n.json ] && return
  echo "$(date -u +%T) start P16 E2 $n" >> $P; cd /root/s
  PYTHONPATH=/root/AS5:/root/AS5/studies/beta_bpeak timeout 7200 python /root/s/fq_check.py --old $OLD --new $NEW --config $c "$@" --json $W/e2/$n.json > $W/e2/$n.log 2>&1
  echo "$(date -u +%T) end P16 E2 $n exit $?" >> $P
}
for q in 2 4 8 16; do fq t3_q$q $C108 --q $q --d_lex 96 --phi 0.55; done
fq n108_q4 $C108 --q 4
fq n108_q8 $C108 --q 8
fq n156_q4 $C156 --q 4
fq n156_q8 $C156 --q 8
fq n108_q16 $C108 --q 16
echo "$(date -u +%T) DONE d23" >> $P
