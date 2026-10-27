#!/bin/bash
# box 11: start stress test 3 (two processes) when both runs of outlier test 3 have ended
P=/root/progress.txt
until grep -q "end K1 H2_d53" $P && grep -q "end K1 H_c63" $P; do sleep 60; done
setsid nohup /root/s/stress3.sh a > /dev/null 2>&1 < /dev/null &
sleep 30
setsid nohup /root/s/stress3.sh b > /dev/null 2>&1 < /dev/null &
