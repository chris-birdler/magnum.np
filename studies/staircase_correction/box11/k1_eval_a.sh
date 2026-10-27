#!/bin/bash
# box 11: K1 evaluation of the A pairs (PLAN 14a T4); the B pairs wait for the cause of the outlier
K="python /root/s/k_eval.py"; R=/root/k1res; mkdir -p $R
A=/root/k1/A; NA="$A/A_n1.pt $A/A_n2.pt $A/A_n3.pt"
pa() { $K --ref $1 --pair $2 $3 $4 --base $5 --ctrl $6 --null $7 $8 $9 --json $R/$2.json > $R/$2.txt 2>&1; echo "$2 exit $?"; }
pa $A/refA.pt A2 $A/A_base.pt $A/A_31a.pt $A/A_base.pt $A/A_ctrl.pt $NA
pa $A/refA.pt A1 $A/A_31a.pt $A/A_a1f.pt $A/A_base.pt $A/A_ctrl.pt $NA
pa $A/refA.pt A3 $A/A_a1f.pt $A/A_3d9.pt $A/A_base.pt $A/A_ctrl.pt $NA
pa $A/refA.pt AB_off $A/A_3d9.pt $A/AB_off.pt $A/A_base.pt $A/A_ctrl.pt $NA
