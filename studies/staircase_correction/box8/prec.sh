#!/bin/bash
# box 8: precision test of PLAN 13 (decision of Chris 2026-10-22, replaces stage 4). Each patch: the S3 run pair (same
# code trees, start state and arguments) again with --precision_cycles double; m after 1, 2, 4, ..., 256 sample
# intervals, and after each interval the tries, accepted steps, t and 2048 probe cells (step_dump.py --trace, as
# /root/prec_step_dump.py). Positive controls: the old codes A_a4 and B_012 with the drive amplitude × (1 + 1e-6).
# A runs (no correction) up to 3 at a time with at least 5 GB free; B runs (correction q 4) alone on the GPU.
# A run that stops with CUDA out of memory starts again (at most 3 tries).
export CUDA_DEVICE=0
P=/root/prof/progress.txt
until [ -f /root/bs3/DONE ]; do sleep 30; done
W=/root/prec; mkdir -p $W
AT="1 2 4 8 16 32 64 128 256"
eval "$(grep '^ARGS=' /root/s3_a1a3f.sh)"; ARGS_A=$ARGS
eval "$(grep '^BARGS=' /root/b_s3.sh)"
I_A=/root/states/K1_2_sat/init.pt
I_B=/root/bs3/prep/init.pt
gpu_free() { nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; }
gpu_apps() { nvidia-smi --query-compute-apps=pid --format=csv,noheader | grep -c .; }
one() { name=$1; code=$2; init=$3; args=$4; D=$W/$name
  for try in 1 2 3; do
    if [ -e $D ] || [ -e $D.log ]; then mkdir -p $W/failed; ts=$(date +%s); mv $D $W/failed/${name}_$ts 2>/dev/null; mv $D.log $W/failed/${name}_$ts.log 2>/dev/null; fi
    mkdir -p $D; echo "$(date -u +%T) start PREC $name try $try" >> $P
    cd /root/$code/studies/beta_bpeak
    PYTHONPATH=/root/$code python /root/prec_step_dump.py --at $AT --trace --save $D/m.pt --go_on -- $args --precision_cycles double --init_from $init --out $D > $D.log 2>&1
    rc=$?
    echo "$(date -u +%T) end PREC $name try $try exit $rc" >> $P
    if [ $rc -eq 0 ] && [ -f $D/m.pt ] && [ -f $D/summary.json ]; then return 0; fi
    grep -q "OutOfMemoryError" $D.log || { echo "$(date -u +%T) FAIL PREC $name (not out of memory)" >> $P; return 1; }
    sleep 120
  done
  echo "$(date -u +%T) FAIL PREC $name (out of memory 3 times)" >> $P; return 1; }
# B5 and B4 run only if their S3 end states in single precision are not bitwise equal
BIT=""; NEED_B5=1; NEED_B4=1
python /root/step_dump.py --compare /root/bs3/b5/old/m.pt /root/bs3/b5/new/m.pt > $W/b5_bitwise.txt 2>&1 && { BIT="$BIT B5"; NEED_B5=0; }
python /root/step_dump.py --compare /root/bs3/b4/old/m.pt /root/bs3/b4/new/m.pt > $W/b4_bitwise.txt 2>&1 && { BIT="$BIT B4"; NEED_B4=0; }
echo "$(date -u +%T) start PREC (bitwise in single precision:$BIT)" >> $P
CTRL_A="${ARGS_A/--b_list 0.006 /--b_list 0.006000006 }"; CTRL_B="${BARGS/--b_list 0.006 /--b_list 0.006000006 }"
[ "$CTRL_A" != "$ARGS_A" ] && [ "$CTRL_B" != "$BARGS" ] || { echo "$(date -u +%T) FAIL PREC control arguments" >> $P; exit 1; }
for x in a4:A_a4 a2:A_a2 a1f:A_a1f a3f:A_a3f ctrlA:A_a4; do
  while [ $(gpu_apps) -ge 3 ] || [ $(gpu_free) -lt 5000 ]; do sleep 15; done
  if [ ${x%%:*} = ctrlA ]; then a="$CTRL_A"; else a="$ARGS_A"; fi
  one ${x%%:*} ${x##*:} $I_A "$a" < /dev/null &
  sleep 90
done
wait
LIST="c63:B_c63 01e:B_01e"
{ [ $NEED_B5 = 1 ] || [ $NEED_B4 = 1 ]; } && LIST="$LIST c56:B_c56"
[ $NEED_B4 = 1 ] && LIST="$LIST d53:B_d53"
LIST="$LIST 012:B_012 df0:B_df0 ctrlB:B_012"
for x in $LIST; do
  while [ $(gpu_apps) -ge 1 ]; do sleep 15; done
  if [ ${x%%:*} = ctrlB ]; then a="$CTRL_B"; else a="$BARGS"; fi
  one ${x%%:*} ${x##*:} $I_B "$a" < /dev/null
done
PAIRS="--pair A2 $W/a4 $W/a2 /root/s3_a2/s3.json $I_A 1e-3 --pair A1 $W/a2 $W/a1f /root/s3_a1f/s3.json $I_A 1e-3 --pair A3 $W/a1f $W/a3f /root/s3_a3f/s3.json $I_A 1e-3 --pair PLAN10 $W/c63 $W/01e /root/bs3/p10.json $I_B 1e-3"
[ $NEED_B5 = 1 ] && PAIRS="$PAIRS --pair B5 $W/01e $W/c56 /root/bs3/b5.json $I_B 1e-3"
[ $NEED_B4 = 1 ] && PAIRS="$PAIRS --pair B4 $W/c56 $W/d53 /root/bs3/b4.json $I_B 1e-3"
PAIRS="$PAIRS --pair B1 $W/012 $W/df0 /root/bs3/b1.json $I_B 0.1"
CTRL="--control A $W/a4 $W/ctrlA 1e-6 /root/s3_a2/s3.json $I_A --control B $W/012 $W/ctrlB 1e-6 /root/bs3/b1.json $I_B"
cd /root; python /root/prec_eval.py $PAIRS $CTRL --bitwise $BIT --json $W/prec.json > $W/prec.txt 2>&1
echo "$(date -u +%T) eval PREC exit $?" >> $P
touch $W/DONE
