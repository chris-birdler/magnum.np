#!/bin/bash
# box 11: after D6, the D7 transfer runs and then the speed measurement (alone on the GPU)
while ! grep -q "DONE d6chain" /root/progress.txt; do sleep 30; done
/root/s/d7t.sh
/root/s/d7.sh
echo "$(date -u +%T) DONE d7chain" >> /root/progress.txt
