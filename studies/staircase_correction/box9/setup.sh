#!/bin/bash
# box 9 (PLAN 14): Python packages, code trees from the bundle, scripts and states. Log: /root/setup.log
set -e
cd /root
pip install -q numpy pyvista scipy setproctitle torchdiffeq xitorch tqdm > /root/pip.log 2>&1
git clone -q /root/up/p14.bundle /root/repo -b main
cd /root/repo
for b in feature/staircase-correction perf/timestep perf/staircase-speed perf/staircase-ba4 perf/speed-all perf/a1-fixed; do
  git fetch -q /root/up/p14.bundle $b:$b
done
for x in 9e0ea7c:M_9e0 9b15ea9:A_9b1 31aa455:A_31a a47492b:A_a1f 3d9de71:A_3d9 7f755a0:AB c63e862:B_c63 01ed90c:B_01e \
         c56c0a5:B_c56 d538bde:B_d53 df0f811:B_df0 b750b1d:BA4; do
  git worktree add -q --detach /root/${x##*:} ${x%%:*}
  echo "${x##*:} $(git -C /root/${x##*:} log --oneline -1 | cut -c1-60)" >> /root/trees.txt
done
mkdir -p /root/s /root/states/K1_2_sat /root/states/K4_2_sat /root/states/Bprep
cp /root/up/scripts/*.py /root/s/
cp /root/up/K1_2_sat_init.pt /root/states/K1_2_sat/init.pt
cp /root/up/K4_2_sat_init.pt /root/states/K4_2_sat/init.pt
cp /root/up/Bprep_init.pt /root/states/Bprep/init.pt
cp /root/up/Bprep_config.json /root/states/Bprep/config.json
cd /root/states && md5sum */init.pt > /root/states.md5
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0)); print('fp64 matmul ok', float((torch.ones(64,64,dtype=torch.float64,device='cuda')@torch.ones(64,64,dtype=torch.float64,device='cuda')).sum()))"
cd /root/AB && CUDA_DEVICE=0 PYTHONPATH=/root/AB python -c "import magnumnp; print('magnumnp import ok')"
echo SETUP DONE
