#!/bin/bash
# PLAN 4.18 step 1 (local CPU): run the named jobs of a job file one after the other, each with OMP_NUM_THREADS
# threads (default 4: two slots on the 8 cores). Output in runs/plan418/<name>, exit codes in runs/plan418/exit_codes.txt.
#   setsid nohup ./run418_slot.sh jobs/plan418_step1.txt B8_s941 timing_T16_s941 > runs/plan418/slot_a.log 2>&1 < /dev/null & disown
JF=$1; shift
cd "$(dirname "$0")" || exit 2
PY=${PY:-python3}
for n in "$@"; do
    line=$(grep "^$n " "$JF") || { echo "$(date '+%F %T') $n: no job line"; exit 2; }
    args=${line#"$n "}
    mkdir -p "runs/plan418/$n"
    echo "$(date '+%F %T') $n start"
    # shellcheck disable=SC2086
    OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} CUDA_DEVICE=-1 "$PY" run_loops.py $args --out "runs/plan418/$n" \
        > "runs/plan418/$n/stdout.log" 2>&1 < /dev/null
    rc=$?
    echo "$(date '+%F %T') $n exit $rc" | tee -a runs/plan418/exit_codes.txt
done
