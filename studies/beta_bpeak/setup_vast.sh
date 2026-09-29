#!/usr/bin/env bash
# Set up a vast.ai instance for the beta(B_peak) study.
#
# Use a PyTorch template image with CUDA torch preinstalled.
# RTX 5090 (Blackwell, sm_120) needs torch built with CUDA >= 12.8.
#
#   curl -sSL https://raw.githubusercontent.com/chris-birdler/magnum.np/study/beta-bpeak-pinning/studies/beta_bpeak/setup_vast.sh | bash
#   or: copy this file to the instance and run  bash setup_vast.sh
set -euo pipefail

REPO=${REPO:-https://github.com/chris-birdler/magnum.np.git}
BRANCH=${BRANCH:-study/beta-bpeak-pinning}
DIR=${DIR:-$HOME/magnum.np}

if [ ! -d "$DIR/.git" ]; then
    git clone --depth 1 -b "$BRANCH" "$REPO" "$DIR"
fi
cd "$DIR"

# do NOT let pip replace the CUDA torch of the image
pip install --no-deps -e .
pip install numpy scipy pyvista setproctitle tqdm torchdiffeq xitorch matplotlib pytest

python - <<'EOF'
import torch, sys
print("torch", torch.__version__, "cuda", torch.version.cuda)
if not torch.cuda.is_available():
    sys.exit("ERROR: no CUDA device visible")
arch = torch.cuda.get_arch_list()
for i in range(torch.cuda.device_count()):
    cap = torch.cuda.get_device_capability(i)
    sm = "sm_%d%d" % cap
    ok = sm in arch or any(a.startswith("compute_") for a in arch)
    print(i, torch.cuda.get_device_name(i), sm, "OK" if ok else "NOT SUPPORTED by this torch build")
    if not ok:
        sys.exit("ERROR: torch build does not support %s (need CUDA >= 12.8 wheels for RTX 50xx)" % sm)
EOF

# quick correctness tests (CPU part + a GPU import)
CUDA_DEVICE=-1 python -m pytest studies/beta_bpeak/test_study.py -q
CUDA_DEVICE=0 python -c "from magnumnp import *; import torch; print('magnum.np on', torch.zeros(1).device)"

python studies/beta_bpeak/jobs.py > /dev/null
echo
echo "Setup complete. Next:"
echo "  cd $DIR/studies/beta_bpeak"
echo "  python run_queue.py jobs/bench.txt --gpus 0"
