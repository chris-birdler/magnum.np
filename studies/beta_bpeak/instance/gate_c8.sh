#!/bin/bash
# Option C, the other 8 corners (Chris 2026-10-18): they start (submit.sh, priority 30) when K3_1_ad has a prepared state
# without a Bloch point (after at most 3 attempts); else runs/ALERT and no start. Writes spool/C8_DECIDED.
#   start: cp gate_c8.sh /root/; setsid nohup bash /root/gate_c8.sh > /root/gate_c8.out 2>&1 < /dev/null & disown
B=/root/magnum.np/studies/beta_bpeak
J=/root/magnum.np-disp/studies/beta_bpeak/jobs
LOG=$B/runs/K3_1_ad/stdout.log
[ -f $B/spool/C8_DECIDED ] && { echo "$(date '+%F %T') already decided"; exit 0; }
until grep -aqE "\[bloch\] attempt [0-9]+: [0-9]+ Bloch points" $LOG 2> /dev/null \
        && { grep -aqE "\[bloch\] attempt [0-9]+: 0 Bloch points|\[bloch\] attempt 2: " $LOG \
             || grep -aq "did not converge" $LOG; } \
        || [ -f $B/runs/K3_1_ad/DONE ]; do
    sleep 120
done
if grep -aqE "\[bloch\] attempt [0-9]+: 0 Bloch points" $LOG; then
    bash /root/submit.sh 30 c8 $J/corners_ad300_rest.txt
    echo "$(date '+%F %T') K3_1_ad clean: $(grep -aE '\[bloch\] attempt' $LOG | tr '\n' ' '); the 8 corners submitted"
else
    echo "$(date '+%F %T') K3_1_ad not clean: $(grep -aE '\[bloch\] attempt' $LOG | tr '\n' ' '); no start" \
        | tee -a $B/runs/ALERT
fi
touch $B/spool/C8_DECIDED
