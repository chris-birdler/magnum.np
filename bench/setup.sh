#!/bin/bash
# Provision a GPU machine (e.g. a vast.ai instance started from a PyTorch
# template image, which ships a CUDA-enabled torch) for benchmarking.
#
# usage: ./setup.sh <git-url-or-local-path> [ref]
set -e

apt-get update -yq && apt-get install -y git
git clone "${1:?usage: ./setup.sh <git-url-or-path> [ref]}" magnum.np
cd magnum.np
git checkout "${2:-main}"
pip install -e .

python - <<'EOF'
import torch
print("torch", torch.__version__, "cuda:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("device:", torch.cuda.get_device_name(0))
EOF
