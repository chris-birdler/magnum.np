#!/bin/bash
# box 11, PLAN 15 (PROTOKOLL 3.44): after the W2 pair FQ (d11.sh FQ ABP ASQ B): (1) the F2 check, a profile of AS2 with
# the source lines (q4_108, 2 profiled sample intervals); (2) W2 for the pairs F1 (ASQ -> AS1, B), F2 (AS1 -> AS2, B
# and A) and F3 (AS2 -> AS3, B and A). A T1 folder is removed when no later pair needs it (disk).
export CUDA_DEVICE=0
P=/root/progress.txt
while ! grep -q "DONE d11 FQ$" $P; do sleep 30; done
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"; B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
W=/root/d10
if [ ! -f $W/f2chk_q4_108.json ]; then
  [ -d /root/kc/c12_AS2_108 ] || cp -a /root/kc/ABP_108 /root/kc/c12_AS2_108
  echo "$(date -u +%T) start D12 f2chk_q4_108" >> $P
  cd /root/AS2/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d12_f2chk PYTHONPATH=/root/AS2 python /root/s/prof_step.py --skip 30 --plain 10 --window 2 --stack --kcache /root/kc/c12_AS2_108 --json $W/f2chk_q4_108.json -- $B --init_from /root/states/Bprep/init.pt --b_list 0.006 --out $W/run_f2chk_q4_108 > $W/f2chk_q4_108.log 2>&1
  rc=$?; echo "$(date -u +%T) end D12 f2chk_q4_108 exit $rc" >> $P
fi
rm -rf /root/w2/t1/B/ABP
bash /root/s/d11.sh F1 ASQ AS1 B; rm -rf /root/w2/t1/B/ASQ
bash /root/s/d11.sh F2 AS1 AS2 B; rm -rf /root/w2/t1/B/AS1
bash /root/s/d11.sh F2A AS1 AS2 A; rm -rf /root/w2/t1/A/AS1
bash /root/s/d11.sh F3 AS2 AS3 B; rm -rf /root/w2/t1/B/AS2
bash /root/s/d11.sh F3A AS2 AS3 A; rm -rf /root/w2/t1/A/AS2
echo "$(date -u +%T) DONE d12" >> $P
