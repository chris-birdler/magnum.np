#!/bin/bash
# box 11: K1 evaluation (PLAN 14a T4). usage: k1_eval.sh pilot | pairs
K="python /root/s/k_eval.py"; R=/root/k1res; mkdir -p $R
A=/root/k1/A; B=/root/k1/B
NA="$A/A_n1.pt $A/A_n2.pt $A/A_n3.pt"; NB="$B/B_n1.pt $B/B_n2.pt $B/B_n3.pt"; NB4="$B/BA4_n1.pt $B/BA4_n2.pt $B/BA4_n3.pt"
if [ "$1" = pilot ]; then
  $K --same $A/refA.pt $A/A_base.pt > $R/same_A.txt 2>&1; echo "same A exit $?"
  $K --same $B/refB.pt $B/B_base.pt > $R/same_B.txt 2>&1; echo "same B exit $? (B is not reproducible on the GPU: information)"
  $K --a4_count $B/refB.pt > $R/a4_count_refB.txt 2>&1; echo "a4 count exit $?"
  $K --ref $A/refA.pt --pilot pilotA --base $A/A_base.pt --ctrl $A/A_ctrl.pt --null $NA --json $R/pilotA.json > $R/pilotA.txt 2>&1; echo "pilot A exit $?"
  $K --ref $B/refB.pt --pilot pilotB --base $B/B_base.pt --ctrl $B/B_ctrl.pt --null $NB --json $R/pilotB.json > $R/pilotB.txt 2>&1; echo "pilot B exit $?"
  $K --ref $B/refB.pt --pilot pilotBA4 --base $B/BA4_base.pt --ctrl $B/BA4_ctrl.pt --null $NB4 --json $R/pilotBA4.json > $R/pilotBA4.txt 2>&1; echo "pilot BA4 exit $?"
else
  pa() { $K --ref $1 --pair $2 $3 $4 --base $5 --ctrl $6 --null $7 $8 $9 ${10} --json $R/$2.json > $R/$2.txt 2>&1; echo "$2 exit $?"; }
  pa $A/refA.pt A2 $A/A_base.pt $A/A_31a.pt $A/A_base.pt $A/A_ctrl.pt $NA
  pa $A/refA.pt A1 $A/A_31a.pt $A/A_a1f.pt $A/A_base.pt $A/A_ctrl.pt $NA
  pa $A/refA.pt A3 $A/A_a1f.pt $A/A_3d9.pt $A/A_base.pt $A/A_ctrl.pt $NA
  pa $A/refA.pt AB_off $A/A_3d9.pt $A/AB_off.pt $A/A_base.pt $A/A_ctrl.pt $NA
  pa $B/refB.pt PLAN10 $B/B_c63.pt $B/B_01e.pt $B/B_base.pt $B/B_ctrl.pt $NB
  pa $B/refB.pt B5 $B/B_01e.pt $B/B_c56.pt $B/B_base.pt $B/B_ctrl.pt $NB
  pa $B/refB.pt B4 $B/B_c56.pt $B/B_base.pt $B/B_base.pt $B/B_ctrl.pt $NB
  pa $B/refB.pt B1 $B/B_base.pt $B/B_df0.pt $B/B_base.pt $B/B_ctrl.pt $NB
  pa $B/refB.pt AB_on $B/BA4_base.pt $B/AB_on.pt $B/BA4_base.pt $B/BA4_ctrl.pt $NB4
  $K --ref $B/refB.pt --pair A4 $B/B_df0.pt $B/BA4_base.pt --base $B/BA4_base.pt --ctrl $B/BA4_ctrl.pt --null $NB4 --align_a4 --json $R/A4.json > $R/A4.txt 2>&1; echo "A4 exit $?"
fi
