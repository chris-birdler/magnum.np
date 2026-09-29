# β(B_peak) study: stress pinning in soft magnetic powder cores

Full micromagnetic AC BH-loops of a periodic FCC powder model (4 spheres,
`DemagFieldPBC`, fp32) to find how particle size, fill factor and
short-range internal stress change the Steinmetz exponent β as a function of
B_peak.

**Start with [HANDOVER.md](HANDOVER.md)** (physics, decisions, procedure,
validation gates, open points).

Quick start on a vast.ai GPU instance:

```bash
bash setup_vast.sh
cd ~/magnum.np/studies/beta_bpeak
python run_queue.py jobs/bench.txt --gpus 0          # timing -> update jobs.py
python run_queue.py jobs/validation.txt --gpus 0,1,2,3
python run_queue.py jobs/production.txt --gpus 0,1,2,3,4,5,6,7
python analyze.py runs/*
```

CPU tests: `CUDA_DEVICE=-1 python -m pytest test_study.py -q`
