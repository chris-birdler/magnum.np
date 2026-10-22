#!/bin/bash
# box 8: S3 of PLAN 11 for the B branch with the correction (q 4) at the geometry of K1_2_sat (N 108):
# PLAN 10 patch (old c63e862, new 01ed90c), B5 (old 01ed90c, new c56c0a5), B4 (old c56c0a5, new d538bde), B1 (old 0125db3, new df0f811)
export CUDA_DEVICE=0
P=/root/prof/progress.txt
until [ -f /root/s3_a3f/DONE ]; do sleep 20; done
W=/root/bs3; mkdir -p $W
for x in c63e862:B_c63 01ed90c:B_01e c56c0a5:B_c56 d538bde:B_d53 0125db3:B_012 df0f811:B_df0; do d=/root/${x##*:}; rm -rf $d; git clone -q /root/perf3.bundle $d -b perf/staircase-speed 2>/dev/null || git clone -q /root/msc8.bundle $d; git -C $d fetch -q /root/msc8.bundle feature/staircase-correction 2>/dev/null; git -C $d checkout -q ${x%%:*}; echo "$d $(git -C $d log --oneline -1 | cut -c1-40)" >> $W/trees.txt; done
BARGS="--d_lex 212.0 --dx_lex 3.0 --phi 0.62 --Leff_lex 30.0 --r_p 1.5 --alpha 0.1 --n_amp 1 --samples 256 --precision single --b_list 0.006 --protocol virgin_asc --init_block 6.0 --relax_mode converge --relax_maxiter 6000 --allow_unconverged --relax_chunk 500 --relax_tol_E 1e-06 --relax_confirm 2 --strayfield newell_ext --strayfield_shells 4 --atol 1e-06 --relax_tol_E_abs 3e-10 --atol_cycles 1e-07 --rtol_cycles 1e-06 --cycles_per_amp 7 --min_calm 2 --max_cycles_per_amp 15 --b_tol 0.01 --seed 2013 --surface fraction --staircase pair --staircase_q 4 --staircase_p 3 --staircase_centroid --max_cycles 1"
# start state: transfer of K1_2_sat onto the fraction surface with the correction, at most 6000 relaxation iterations
echo "$(date -u +%T) start B S3 prep" >> $P
cd /root/B_df0/studies/beta_bpeak; PYTHONPATH=/root/B_df0 python /root/step_dump.py --at 1 --save $W/prep_m.pt -- $BARGS --init_transfer /root/states/K1_2_sat/init.pt --out $W/prep > $W/prep.log 2>&1
echo "$(date -u +%T) end B S3 prep exit $?" >> $P
I0=$W/prep/init.pt
for k in 1 2 3; do PYTHONPATH=/root/B_df0 python /root/perturb_init.py $I0 $W/null$k/init.pt --seed $k > $W/perturb_$k.txt 2>&1; done
jobs=$W/jobs.txt; : > $jobs
for set in "p10 B_c63 B_01e" "b5 B_01e B_c56" "b4 B_c56 B_d53" "b1 B_012 B_df0"; do
  set -- $set
  echo "$1 old $2 $I0" >> $jobs
  for k in 1 2 3; do echo "$1 null$k $2 $W/null$k/init.pt" >> $jobs; done
  echo "$1 new $3 $I0" >> $jobs
done
one() { s=$1; name=$2; code=$3; init=$4; D=$W/$s/$name; mkdir -p $D; echo "$(date -u +%T) start BS3 $s $name" >> $P; cd /root/$code/studies/beta_bpeak; PYTHONPATH=/root/$code python /root/step_dump.py --at 256 --save $D/m.pt --go_on -- $BARGS --init_from $init --out $D > $D.log 2>&1; echo "$(date -u +%T) end BS3 $s $name exit $?" >> $P; }
export -f one; export W P BARGS CUDA_DEVICE
i=0; while read s name code init; do
  while [ $(jobs -rp | wc -l) -ge 3 ]; do sleep 10; done
  one $s $name $code $init &
  i=$((i+1)); sleep 150                                   # the set-ups of the correction (up to about 6 GB) do not overlap
done < $jobs
wait
for s in p10 b5 b4 b1; do python /root/s3_eval.py --old $W/$s/old --null $W/$s/null1 $W/$s/null2 $W/$s/null3 --new $W/$s/new --init_old $I0 --init_null $W/null1/init.pt $W/null2/init.pt $W/null3/init.pt --json $W/$s.json > $W/$s.txt 2>&1; done
touch $W/DONE
