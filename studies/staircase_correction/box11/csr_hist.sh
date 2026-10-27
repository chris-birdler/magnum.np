#!/bin/bash
# box 11: row lengths of the CSR matrix of the correction (geometry B, d538bde) on the CPU, beside the GPU runs; the event indices of PLAN 14a
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
mkdir -p /root/csr
echo "$(date -u +%T) start csr_hist" >> $P; cd /root/B_d53/studies/beta_bpeak
CUDA_DEVICE=-1 OMP_NUM_THREADS=16 TORCHINDUCTOR_CACHE_DIR=/root/ic/csr PYTHONPATH=/root/B_d53 timeout 5400 nice -n 19 python /root/s/csr_hist.py 1513981 876667 1634169 3060786 3615150 3632484 2605404 -- $ARGS_B --b_list 0.006 --out /root/csr/run > /root/csr/csr_hist.log 2>&1
rc=$?; echo "$(date -u +%T) end csr_hist exit $rc" >> $P
