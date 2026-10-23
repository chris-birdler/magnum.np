#!/bin/bash
# box 11, phase 1 of PLAN 14a: K1 references, then the pilot arms (2 runs at a time on the GPU)
/root/s/k1.sh /root/s/p1_refs.txt 2
/root/s/k1.sh /root/s/p1_pilot.txt 2
touch /root/k1/P1_DONE
