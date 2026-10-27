#!/bin/bash
# box 11: after the D5 arms, the pair evaluation of D5 (CPU) and D6 (t12p.sh, GPU) at the same time
while ! grep -q "DONE /root/s/p3_arms.txt" /root/progress.txt; do sleep 30; done
/root/s/k1_eval_p.sh pairs > /root/k1res_p/pairs.out 2>&1 &
/root/s/t12p.sh
wait
echo "$(date -u +%T) DONE d6chain" >> /root/progress.txt
