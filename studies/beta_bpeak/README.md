# β(B_peak) study: what controls the Steinmetz exponent in powder cores

Full micromagnetic AC BH-loops of a periodic FCC powder model (4 spheres,
`DemagFieldPBC`, fp32) for a fictitious material in reduced units. Scan of
pinning strength κ, particle size d/l_ex, fill factor φ and anisotropy Q.

**Start with [HANDOVER.md](HANDOVER.md).**

```bash
bash setup_vast.sh
cd ~/magnum.np/studies/beta_bpeak
python run_queue.py jobs/bench.txt --gpus 0          # timing -> update jobs.py
python run_queue.py jobs/production.txt --gpus 0,1,2,3
python run_queue.py jobs/control.txt --gpus 4,5,6,7
python analyze.py runs/*
```

CPU tests: `CUDA_DEVICE=-1 python -m pytest test_study.py -q`
