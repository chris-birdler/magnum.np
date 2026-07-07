# magnum.np GPU benchmarks

Self-contained benchmark scripts to compare two versions of magnum.np on real
hardware (e.g. a rented vast.ai GPU instance). The scripts only depend on an
installed `magnumnp` and the standard library, so they can benchmark any
checked-out ref.

## Procedure (vast.ai)

1. Rent an instance (e.g. RTX 4090) from a **PyTorch template image** (CUDA
   torch preinstalled) and copy this repository (or just clone it) onto the
   machine:

   ```bash
   ./bench/setup.sh https://gitlab.com/magnum.np/magnum.np.git <feature-branch>
   cp -r magnum.np/bench ~/bench   # keep the scripts checkout-independent
   ```

2. Baseline run:

   ```bash
   cd magnum.np && git checkout <baseline-sha> && pip install -e . && cd ..
   ~/bench/run_all.sh results/base
   ```

3. Optimized run:

   ```bash
   cd magnum.np && git checkout <feature-branch> && cd ..
   ~/bench/run_all.sh results/opt
   ```

4. Compare (exits non-zero if anything regresses by more than 3%):

   ```bash
   python ~/bench/compare.py results/base results/opt
   ```

5. Run the unit tests once on the GPU (the CI has no GPU runner):

   ```bash
   cd magnum.np && python -m pytest tests/unit
   ```

## What the benchmarks measure

| script               | metric                          | what it validates                                    |
|----------------------|---------------------------------|------------------------------------------------------|
| `bench_kernel_init`  | init wall time, GPU/CPU peak    | chunked CPU-resident kernel setup (expect the GPU peak to drop from 10-20x kernel size to ~the kernel size) |
| `bench_demag_h`      | ms per `DemagField.h()`         | per-component FFT loop vs batched single-FFT (adopt batched only if >=10% faster at *both* precisions) |
| `bench_sp4`          | LLG steps/s, peak memory        | end-to-end effect of the solver/field-term optimizations; f32 vs f64 speedup |
| `bench_minimize`     | MinimizerBB iterations/s        | compiled `_midpoint`, clone removal, sync reduction   |

Every script runs each variant in a fresh subprocess (clean allocator/peak
statistics, isolated dtype) and writes one JSON file per benchmark into the
given results directory. All scripts also run on CPU-only machines (reduced
iteration counts, `CUDA_DEVICE=-1`).

## Acceptance criteria for the optimization branch

- `kernel_init` GPU peak reduced by >= 5x
- `sp4` f64 steps/s at least at baseline (target +10..30%)
- `sp4` f32 >= 1.5x faster than f64
- no benchmark regresses by more than 3%
