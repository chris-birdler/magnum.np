#!/bin/bash
# box 11: after the D5 pilot, the pilot evaluation (CPU) and check X1 (GPU) at the same time
while ! grep -q "DONE /root/s/p3_pilot.txt" /root/progress.txt; do sleep 30; done
mkdir -p /root/k1res_p
/root/s/k1_eval_p.sh pilot > /root/k1res_p/pilot.out 2>&1 &
/root/s/x1.sh
wait
echo "$(date -u +%T) DONE x1chain" >> /root/progress.txt
