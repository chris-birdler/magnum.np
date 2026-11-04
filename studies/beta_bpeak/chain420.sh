#!/bin/bash
# PLAN 4.20 box chain (go of Chris 2026-11-04): stage B, the selection, stage C, on the dispatcher of the box. Runs on
# the box without the laptop: setsid nohup bash chain420.sh <gate A variants, e.g. V1,V2> <q16 mode> > chain420.out ...
# q16 mode: "serial" (host memory below 128 GB: one q-16 job at a time) or "setup" (128 GB or more: the next q-16 job
# starts when the set-up of the one before has ended). Log: runs/chain420.log; end marker runs/CHAIN420_END (text: the
# reason). Gates of PLAN 4.20: stage B -> eval_420.py stageB -> SELECTED (V1, V2 or NONE) -> stage C of that variant;
# NONE (or no selection): runs/STOP_NEW_JOBS at once, no stage C (gate 5: destroy). After stage C: eval_420.py stageC;
# with fewer than 7 counted q-8 jobs that ran their cycles, the spare seed 925 runs and the evaluation is made again.
# Control R0 (a background check): a fail writes runs/R0_FAIL, stops new jobs and ends the chain (the jobs stop).
# Priorities (smaller first): stage B 10, tail of stage B 11 (R0, HFP_s922), q-16 jobs 12 (one at a time), stage C 15,
# spare 25.
GATE=$1
MODE=$2
B=${CHAIN_B:-/root/magnum.np/studies/beta_bpeak}     # CHAIN_B, CHAIN_SUBMIT, CHAIN_POLL: only for the local test
SUBMIT=${CHAIN_SUBMIT:-/root/submit.sh}
POLL=${CHAIN_POLL:-120}
J=$B/jobs
S=$B/spool
R=$B/runs
L=$R/chain420.log
cd $B || exit 1
[ -n "$GATE" ] && [ "$MODE" = serial -o "$MODE" = setup ] || { echo "usage: chain420.sh <V1,V2> <serial|setup>"; exit 1; }
log() { echo "$(date -u '+%F %T') $*" >> $L; }
stop_chain() {   # stop_chain <reason>: no new job starts (the dispatcher's stop file), end marker, end of the chain
    touch $R/STOP_NEW_JOBS
    echo "$1" > $R/CHAIN420_END
    log "STOP: $1"
    kill $(jobs -p) 2> /dev/null
    exit 0
}

submit() {   # submit <priority> <tag> <job file> [names ...]: only the named lines (all without names)
    local p=$1 t=$2 f=$3; shift 3
    local tmp=$R/.sub_$t
    if [ $# -eq 0 ]; then grep -v '^#' "$f" | grep . > $tmp; else
        : > $tmp; for n in "$@"; do grep "^$n " "$f" >> $tmp; done; fi
    bash $SUBMIT "$p" "$t" "$tmp" >> $L 2>&1
}
wait_done() {   # wait_done <priority> <tag>: until the dispatcher wrote the .done file of the spool file
    while [ ! -f "$S/$1_$2.jobs.done" ]; do
        [ -f $R/R0_FAIL ] && stop_chain "R0 failed: $(cat $R/R0_FAIL)"
        sleep $POLL; done
    log "done $1_$2: $(cat $S/$1_$2.jobs.done | head -c 300)"
}
q16_chain() {   # q16_chain <tag> <job file>: the q-16 jobs one by one (serial) or with one set-up at a time
    local t=$1 f=$2 n
    for n in $(grep -v '^#' "$f" | cut -d' ' -f1); do
        submit 12 "${t}_$n" "$f" "$n"
        if [ "$MODE" = serial ]; then
            wait_done 12 "${t}_$n"
        else
            while [ ! -f "$S/12_${t}_$n.jobs.done" ] && ! grep -q "^\[relax\] it\|^\[relax_try" "$R/$n/stdout.log" 2>/dev/null; do
                [ -f $R/R0_FAIL ] && stop_chain "R0 failed: $(cat $R/R0_FAIL)"
                sleep $POLL; done
            log "q16 $n: set-up ended"
        fi
    done
}

log "chain420 start: gate A variants $GATE, q16 mode $MODE"
# stage B: F8 first, then D8, the variants of gate A, C8
names=""
for a in F8 D8 $(echo "$GATE" | tr ',' ' ' | sed 's/V/D8V/g') C8; do
    for s in 941 942 943 944; do names="$names ${a}_s$s"; done
done
submit 10 stageB $J/plan420_stageB.txt $names
# the tail of stage B (independent of the selection): R0, HFP_s922, then the D16 jobs
submit 11 tailB $J/plan420_tailB.txt
( while [ ! -f $R/R0_s921/DONE ] && [ ! -f $S/11_tailB.jobs.done ]; do sleep $POLL; done
  python eval_420.py r0 $R --ref416 $B/results/hf416.json >> $L 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then echo "eval_420.py r0 exit $rc (see chain420.log)" > $R/R0_FAIL; fi ) &
q16_chain d16 $J/plan420_q16_D16.txt &
QP=$!
wait_done 10 stageB
python eval_420.py stageB $R --gateA "$GATE" --out $R/plan420_stageB >> $L 2>&1
SEL=$(cat $R/SELECTED 2>/dev/null)
log "selection: ${SEL:-none (eval failed)}"
if [ "$SEL" != V1 ] && [ "$SEL" != V2 ]; then
    stop_chain "no stage C: selection ${SEL:-eval failed}"
fi
submit 15 stageC_$SEL $J/plan420_stageC_$SEL.txt
wait $QP                                    # the D16 chain before the Dv16 chain (one q-16 set-up at a time)
q16_chain dv16 $J/plan420_q16_$SEL.txt
wait_done 15 stageC_$SEL
wait_done 11 tailB
for n in $(grep -v '^#' $J/plan420_q16_D16.txt | cut -d' ' -f1); do wait_done 12 "d16_$n"; done
for n in $(grep -v '^#' $J/plan420_q16_$SEL.txt | cut -d' ' -f1); do wait_done 12 "dv16_$n"; done
[ -f $R/R0_FAIL ] && stop_chain "R0 failed: $(cat $R/R0_FAIL)"
python eval_420.py stageC $R --variant $SEL --ref416 $B/results/hf416.json --out $R/plan420_stageC >> $L 2>&1
if python - $R/plan420_stageC.json <<'EOF'
import json, sys
ca = json.load(open(sys.argv[1]))["Ca"]["q8"]
sys.exit(0 if ca["status"] == "not determined" and len(ca["ran"]) < 7 else 1)
EOF
then
    log "fewer than 7 counted q-8 jobs ran their cycles: the spare seed 925"
    submit 25 spare_$SEL $J/plan420_spare_$SEL.txt
    wait_done 25 spare_$SEL
    python eval_420.py stageC $R --variant $SEL --ref416 $B/results/hf416.json --out $R/plan420_stageC >> $L 2>&1
fi
echo "stage C done: variant $SEL" > $R/CHAIN420_END
log "end"
