#!/bin/bash
# box 11: start the second process of stress test 3 when H_c63 has ended (the first process started at the end of H2_d53)
P=/root/progress.txt
until grep -q "end K1 H_c63" $P; do sleep 60; done
setsid nohup /root/s/stress3.sh b > /dev/null 2>&1 < /dev/null &
