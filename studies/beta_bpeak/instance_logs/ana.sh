#!/bin/bash
# domain analysis of the snapshot runs, in its own checkout (the chain keeps its code)
set -e
cd /root
[ -d ana ] || git clone -q --depth 1 -b study/beta-bpeak-pinning https://github.com/chris-birdler/magnum.np.git ana
cd /root/ana && git fetch -q --depth 1 origin study/beta-bpeak-pinning && git checkout -q -f FETCH_HEAD
cd /root/ana/studies/beta_bpeak
R=/root/magnum.np/studies/beta_bpeak/runs
git log --oneline -1 > /root/ana.log
for r in V_s901 V_s902 B300_s921 B300_s922 B300_s923 B300_s924 B300_s925 B300_s926 B300_s927 B300_s928; do
  while [ ! -f $R/$r/DONE ]; do sleep 60; done
  echo $R/$r
done | xargs -P 10 -I{} env OMP_NUM_THREADS=3 PYTHONPATH=/root/ana python domains_batch.py {} >> /root/ana.log 2>&1
touch /root/ANA_DONE
