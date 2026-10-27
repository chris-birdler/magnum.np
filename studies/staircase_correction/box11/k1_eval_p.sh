#!/bin/bash
# box 11: K1 evaluation of the B side with the P trees (PLAN 14b D5). usage: k1_eval_p.sh pilot | pairs
K="python /root/s/k_eval.py"; R=/root/k1res_p; mkdir -p $R
B=/root/k1/B
NB="$B/BP_n1.pt $B/BP_n2.pt $B/BP_n3.pt"; NB4="$B/BA4P_n1.pt $B/BA4P_n2.pt $B/BA4P_n3.pt"
if [ "$1" = pilot ]; then
  $K --same $B/BP_base.pt $B/BP_base2.pt > $R/same_BP.txt 2>&1; echo "same BP base/base2 exit $?"
  $K --ref $B/refB.pt --pilot pilotBP --base $B/BP_base.pt --ctrl $B/BP_ctrl.pt --null $NB --json $R/pilotBP.json > $R/pilotBP.txt 2>&1; echo "pilot BP exit $?"
  $K --ref $B/refB.pt --pilot pilotBA4P --base $B/BA4P_base.pt --ctrl $B/BA4P_ctrl.pt --null $NB4 --json $R/pilotBA4P.json > $R/pilotBA4P.txt 2>&1; echo "pilot BA4P exit $?"
else
  pa() { $K --ref $1 --pair $2 $3 $4 --base $5 --ctrl $6 --null $7 $8 $9 ${10} --json $R/$2.json > $R/$2.txt 2>&1; echo "$2 exit $?"; }
  pa $B/refB.pt PLAN10 $B/B_c63.pt $B/B_01eP.pt $B/BP_base.pt $B/BP_ctrl.pt $NB
  pa $B/refB.pt B5 $B/B_01eP.pt $B/B_c56P.pt $B/BP_base.pt $B/BP_ctrl.pt $NB
  pa $B/refB.pt B4 $B/B_c56P.pt $B/BP_base.pt $B/BP_base.pt $B/BP_ctrl.pt $NB
  pa $B/refB.pt B1 $B/BP_base.pt $B/B_df0P.pt $B/BP_base.pt $B/BP_ctrl.pt $NB
  pa $B/refB.pt AB_on $B/BA4P_base.pt $B/ABP_on.pt $B/BA4P_base.pt $B/BA4P_ctrl.pt $NB4
  $K --ref $B/refB.pt --pair A4 $B/B_df0P.pt $B/BA4P_base.pt --base $B/BA4P_base.pt --ctrl $B/BA4P_ctrl.pt --null $NB4 --align_a4 --json $R/A4.json > $R/A4.txt 2>&1; echo "A4 exit $?"
fi
