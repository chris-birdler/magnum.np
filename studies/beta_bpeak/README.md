# β(B_peak) study: what controls the Steinmetz exponent in powder cores

Full micromagnetic AC BH-loops of a periodic FCC powder model (4 spheres of
1 µm, `DemagFieldPBC`, fp32) for a fictitious nanocrystalline material in
reduced units, for B_peak = 10 … 150 mT. The anisotropy is the Herzer residual anisotropy (cubes of
edge L_eff with random axes, Q_eff = (l_ex/L_eff)²) plus an optional residual
stress anisotropy per particle (r_p = K_p/K_eff). Factors: Q_eff, r_p, φ.

**Start with [PROTOKOLL.md](PROTOKOLL.md).**

```bash
bash setup_vast.sh
cd ~/magnum.np/studies/beta_bpeak
python run_queue.py jobs/bench.txt --gpus 0          # timing -> update jobs.py
python run_queue.py jobs/meshtest.txt --gpus 0,1,2,3 # FIRST: grid pinning of vortex cores? (PROTOKOLL 7.1)
python analyze.py --compare runs/M_L12_dx3 runs/M_L12_dx15
python run_queue.py jobs/steady.txt --gpus 0,1       # steady cycles (PROTOKOLL 7.2)
python run_queue.py jobs/numfloor.txt --gpus 0,1     # numerical floor of small loops (PROTOKOLL 7.3)
python run_queue.py jobs/pilot.txt --gpus 0,1,2,3    # pilot, then Chris decides the production matrix
python analyze.py runs/T_*
python analyze.py --pair runs/T_L18 runs/T_L18_f2
```

CPU tests: `CUDA_DEVICE=-1 python -m pytest test_study.py -q`
