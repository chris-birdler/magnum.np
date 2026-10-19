#!/bin/bash
# usage: ./run_all.sh results/<label>
set -e
OUT=${1:?usage: ./run_all.sh results/<label>}
DIR=$(dirname "$0")

python "$DIR/bench_kernel_init.py" --out "$OUT"
python "$DIR/bench_demag_h.py"     --out "$OUT"
python "$DIR/bench_sp4.py"         --out "$OUT"
python "$DIR/bench_minimize.py"    --out "$OUT"

echo "results written to $OUT"
