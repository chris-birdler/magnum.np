#!/bin/bash
# box 11, PLAN 14a T4 (K1): one job line = NAME TREE GEO B_LIST STEP_DUMP_OPTIONS...; GEO A or B; up to $PAR at a time.
# usage: k1.sh JOBFILE PAR
export CUDA_DEVICE=0
P=/root/progress.txt
ARGS_A="--d_lex 212.0 --dx_lex 3.0 --phi 0.62 --Leff_lex 30.0 --r_p 1.5 --alpha 0.1 --n_amp 1 --samples 256 --precision single --protocol virgin_asc --init_block 6.0 --relax_mode converge --relax_maxiter 80000 --relax_chunk 500 --relax_tol_E 1e-06 --relax_confirm 2 --strayfield newell_ext --strayfield_shells 4 --atol 1e-06 --relax_tol_E_abs 3e-10 --atol_cycles 1e-07 --rtol_cycles 1e-06 --cycles_per_amp 7 --min_calm 2 --max_cycles_per_amp 15 --b_tol 0.01 --init_mode saturated --predemag_h 0 --h_sat 2 --bloch_free_tries 3 --seed 2013 --sat_static --max_cycles 1 --precision_cycles double --init_from /root/states/K1_2_sat/init.pt"
ARGS_B="--d_lex 212.0 --dx_lex 3.0 --phi 0.62 --Leff_lex 30.0 --r_p 1.5 --alpha 0.1 --n_amp 1 --samples 256 --precision single --protocol virgin_asc --init_block 6.0 --relax_mode converge --relax_maxiter 6000 --allow_unconverged --relax_chunk 500 --relax_tol_E 1e-06 --relax_confirm 2 --strayfield newell_ext --strayfield_shells 4 --atol 1e-06 --relax_tol_E_abs 3e-10 --atol_cycles 1e-07 --rtol_cycles 1e-06 --cycles_per_amp 7 --min_calm 2 --max_cycles_per_amp 15 --b_tol 0.01 --seed 2013 --surface fraction --staircase pair --staircase_q 4 --staircase_p 3 --staircase_centroid --max_cycles 1 --precision_cycles double --init_from /root/states/Bprep/init.pt"
one() {
  name=$1; tree=$2; geo=$3; bl=$4; shift 4
  W=/root/k1/$geo; mkdir -p $W
  [ -f $W/$name.pt ] && { echo "$(date -u +%T) skip K1 $name (done)" >> $P; return; }
  echo "$(date -u +%T) start K1 $name" >> $P
  cd /root/$tree/studies/beta_bpeak
  if [ $geo = A ]; then A="$ARGS_A"; else A="$ARGS_B"; fi
  TORCHINDUCTOR_CACHE_DIR=/root/ic/$geo/$name PYTHONPATH=/root/$tree python /root/s/step_dump.py "$@" -- $A --b_list $bl --out $W/run_$name > $W/$name.log 2>&1
  rc=$?; echo "$(date -u +%T) end K1 $name exit $rc" >> $P
}
export -f one; export P ARGS_A ARGS_B CUDA_DEVICE
grep -v '^#' $1 | grep -v '^$' | xargs -P ${2:-3} -L1 bash -c 'one "$@"' _
echo "$(date -u +%T) DONE $1" >> $P
