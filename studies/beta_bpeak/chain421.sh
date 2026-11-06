#!/bin/bash
# PLAN 4.21 box chain (gap test of V2 with the production line): the 12 runs of jobs/plan421.txt on the dispatcher of
# the box, with the checks of the plan. Runs on the box without the laptop:
#   setsid nohup bash chain421.sh <q16 mode> <d5> > /root/chain421.out 2>&1 < /dev/null & disown
# q16 mode: "setup" (host memory 128 GB or more: the next job is submitted when the set-up of the one before has
# ended, so that one q-16 set-up of about 34.5 GiB runs at a time) or "serial" (one job at a time).
# d5 (decision D5 of Chris): "a" = a Bloch sample in a run of K1, K3 or K4 stops new jobs; "b" = all jobs continue.
# Order (PLAN 4.21): the small runs K1_0, K3_0, K1_1, K3_1 first; the large runs only after the first-run check (R4) of
# the first finished run passed; then the rest in the order of the job file. If all submitted small runs end and none
# passes R4, the chain stops new jobs.
# For each finished run (DONE), in the background: check_sat_run.py (R4; its last line "C5 data part: PASS|FAIL");
# after a PASS: eval_hotspots.py, eval_hotcores.py and census_bloch.py; only if all three end with exit code 0,
# instance_postproc.py --once for this run (it deletes the snapshots and the accumulator files). Else the files stay
# (fetch421.sh --kept downloads them). A run with R4_RESULT is not checked again (a restart of the chain).
# Log: runs/chain421.log; end marker runs/CHAIN421_END (text: the reason); a stop: runs/STOP_NEW_JOBS.
MODE=$1
D5=$2
B=${CHAIN_B:-/root/magnum.np/studies/beta_bpeak}     # CHAIN_B, CHAIN_SUBMIT, CHAIN_POLL, CHAIN_PY: for the local test
SUBMIT=${CHAIN_SUBMIT:-/root/submit.sh}
POLL=${CHAIN_POLL:-120}
PY=${CHAIN_PY:-python}
J=$B/jobs/plan421.txt
S=$B/spool
R=$B/runs
L=$R/chain421.log
SMALL="K1_0_v16 K3_0_v16 K1_1_v16 K3_1_v16"
SETUP_END='^\[relax\] virgin state\|^\[relax\] it\|^\[relax_try'   # first lines after the set-up of run_loops
cd $B || exit 1
[ "$MODE" = serial -o "$MODE" = setup ] && [ "$D5" = a -o "$D5" = b ] || {
    echo "usage: chain421.sh <serial|setup> <a|b>"; exit 1; }
mkdir -p $R
log() { echo "$(date -u '+%F %T') $*" >> $L; }
NAMES=$(grep -v '^#' $J | grep . | cut -d' ' -f1)

stopped() { [ -f $R/STOP_NEW_JOBS ]; }
stop_new() {   # stop_new <reason>: no new job starts; the jobs in progress end normally
    stopped && return
    echo "$1" > $R/STOP_REASON
    touch $R/STOP_NEW_JOBS
    log "STOP NEW JOBS: $1"
}

submit_one() {   # submit_one <name>: the line of the job file as its own spool file (priority 10)
    local n=$1 tmp=$R/.sub_$1
    grep "^$n " $J > $tmp
    bash $SUBMIT 10 "$n" "$tmp" >> $L 2>&1
}

done_jobs() { ls $S/10_*.jobs.done 2> /dev/null | sed 's#.*/10_##; s#\.jobs\.done##'; }
submitted() { ls $S/10_*.jobs 2> /dev/null | sed 's#.*/10_##; s#\.jobs##'; }
r4_pass_first() { for n in $(submitted); do [ "$(cat $R/$n/R4_RESULT 2> /dev/null)" = PASS ] && { echo $n; return; }; done; }

bloch_k134() {   # the runs of K1, K3 or K4 with a Bloch sample in a cycle (bloch_<i>_amp<i>.csv: columns 3 and 4)
    local f
    for f in $R/K[134]_*_v16/bloch_[0-9][0-9]_amp[0-9][0-9].csv; do
        [ -f "$f" ] || continue
        awk -F, '!/^#/ && ($3 + $4) > 0 {found = 1} END {exit !found}' "$f" && echo "$(basename $(dirname $f))"
    done | sort -u
}

postproc_run() {   # postproc_run <name>: R4, then (PASS) the evaluations and the post-processing of one run
    local n=$1 res rc1 rc2 rc3
    res=$($PY check_sat_run.py $R/$n > $R/$n/check_sat_run.txt 2>&1; tail -1 $R/$n/check_sat_run.txt)
    log "$n R4: $res ($(grep -c '^FAIL' $R/$n/check_sat_run.txt) FAIL lines)"
    if [ "$res" != "C5 data part: PASS" ]; then
        echo FAIL > $R/$n/R4_RESULT
        stop_new "R4 FAIL in $n"
        return
    fi
    echo PASS > $R/$n/R4_RESULT
    $PY eval_hotspots.py --json $R/$n/hotspots.json $R/$n > $R/$n/hotspots.md 2> $R/$n/hotspots.err; rc1=$?
    $PY eval_hotcores.py --json $R/$n/hotcores.json $R/$n > $R/$n/hotcores.md 2> $R/$n/hotcores.err; rc2=$?
    $PY census_bloch.py --out $R/$n/census_bloch.json $R/$n > $R/$n/census_bloch.log 2>&1; rc3=$?
    log "$n eval_hotspots exit $rc1, eval_hotcores exit $rc2, census_bloch exit $rc3"
    if [ $rc1 -ne 0 ] || [ $rc2 -ne 0 ] || [ $rc3 -ne 0 ]; then
        log "$n: an evaluation failed; the snapshots and the accumulator files stay (fetch421.sh --kept)"
        echo "evaluation exit $rc1 $rc2 $rc3" > $R/$n/POSTPROC_KEPT
        return
    fi
    $PY instance_postproc.py --runs $R --stop $R/.no_stop --once --workers 1 --pattern "^$n\$" \
        >> $R/postproc421.log 2>&1
    log "$n post-processing exit $?, POSTPROC_DONE $([ -f $R/$n/POSTPROC_DONE ] && echo yes || echo NO)"
}

declare -A HANDLED
handle_finished() {   # start postproc_run in the background for each run that ended since the last call
    local n
    for n in $(done_jobs); do
        [ -n "${HANDLED[$n]}" ] && continue
        HANDLED[$n]=1
        if [ -f $R/$n/R4_RESULT ]; then
            log "$n: R4_RESULT exists (checked before a restart of the chain)"
        elif [ ! -f $R/$n/DONE ]; then
            log "$n ended without DONE (failed job): $(tail -c 300 $S/10_$n.jobs.done | tr '\n' ' ')"
        else
            postproc_run $n &
        fi
    done
}

checks() {   # all checks of one poll
    handle_finished
    if [ "$D5" = a ]; then
        local b
        b=$(bloch_k134 | tr '\n' ' ')
        [ -n "$b" ] && stop_new "D5 (a): a Bloch sample in $b"
    fi
}

wait_setup() {   # wait_setup <name>: until its set-up has ended (or the job ended); serial: until the job ended.
    local n=$1   # A job that did not start before a stop never starts: return then.
    while [ ! -f "$S/10_$n.jobs.done" ]; do
        if [ "$MODE" = setup ] && grep -q "$SETUP_END" "$R/$n/stdout.log" 2> /dev/null; then
            log "$n: set-up ended"
            return
        fi
        if stopped && [ ! -f "$R/$n/stdout.log" ]; then
            log "$n: not started before the stop"
            return
        fi
        checks
        sleep $POLL
    done
}

wait_r4() {   # until a run passed R4, or a stop, or all submitted jobs ended without an R4 PASS
    local all
    while [ -z "$(r4_pass_first)" ] && ! stopped; do
        checks
        all=1
        for n in $(submitted); do [ -f $S/10_$n.jobs.done ] && [ -f $R/$n/R4_RESULT -o ! -f $R/$n/DONE ] || all=0; done
        if [ $all = 1 ] && [ -z "$(r4_pass_first)" ]; then
            stop_new "no small run passed R4 (all submitted runs ended)"
        fi
        sleep $POLL
    done
}

log "chain421 start: q16 mode $MODE, D5 $D5, $(echo $NAMES | wc -w) jobs"
for n in $NAMES; do
    stopped && break
    if ! echo " $SMALL " | grep -q " $n "; then        # a large run (or a later small run): only after R4 PASS
        wait_r4
        stopped && break
    fi
    submit_one $n
    log "submitted $n"
    wait_setup $n
done
# the jobs in progress end (the jobs that did not start before a stop stay in spool/); then the last checks
while true; do
    checks
    running=$(for n in $(submitted); do [ -f $S/10_$n.jobs.done ] && continue
                  if ! stopped || [ -f $R/$n/stdout.log ]; then echo $n; fi; done | wc -l)
    [ "$running" -eq 0 ] && break
    sleep $POLL
done
checks
wait                                           # the post-processing in the background
echo "$(stopped && echo "stopped: $(cat $R/STOP_REASON); ")chain421 end: $(done_jobs | wc -l) jobs ended" > $R/CHAIN421_END
log "end: $(cat $R/CHAIN421_END)"
