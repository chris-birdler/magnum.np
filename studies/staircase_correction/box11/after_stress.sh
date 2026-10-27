#!/bin/bash
# box 11: start the A patch arms when the stress test is done
until grep -q "DONE stress" /root/progress.txt; do sleep 30; done
/root/s/k1.sh /root/s/p2a_arms.txt 2
