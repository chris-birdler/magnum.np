#!/bin/bash
# box 11, PLAN 15 WQ (PROTOKOLL 3.44), after d17.sh: diagnosis of the stop of the transfer T16_q16_156 (exit code 137,
# the OOM killer of the container, limit 89.0e9 bytes): mem_phases.py gives the host memory in each phase of the geometry
# and of the set-up at q 16, N 156 (tree AS3). The cgroup counters before and after.
export CUDA_DEVICE=0
P=/root/progress.txt
while ! grep -q "DONE d17$" $P; do sleep 30; done
Q=/root/wq; C=/sys/fs/cgroup/memory
echo "before: oom_kill $(grep oom_kill $C/memory.oom_control | tail -1), max_usage $(cat $C/memory.max_usage_in_bytes)" > $Q/mem_q16_156_cgroup.txt
echo "$(date -u +%T) start D18 mem q16_156" >> $P
cd /root/s
PYTHONPATH=/root/AS3:/root/AS3/studies/beta_bpeak timeout 3600 python /root/s/mem_phases.py --config /root/states/T16_q16_156/config.json > $Q/mem_q16_156.log 2>&1
rc=$?; echo "$(date -u +%T) end D18 mem q16_156 exit $rc" >> $P
echo "after: oom_kill $(grep oom_kill $C/memory.oom_control | tail -1), max_usage $(cat $C/memory.max_usage_in_bytes)" >> $Q/mem_q16_156_cgroup.txt
echo "$(date -u +%T) DONE d18" >> $P
