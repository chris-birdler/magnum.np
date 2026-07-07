"""End-to-end LLG stepping throughput on sp4-like problems (f32/f64).

Captures the combined effect of the solver/field-term optimizations and
validates the single-precision speedup on real hardware.
"""
import json
import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SIZES = [(200, 50, 1), (500, 125, 2)]


def run_variant(dtype_arg, nx, ny, nz):
    import torch
    import magnumnp as mag
    from common import reset_peak_memory, peak_memory, set_dtype

    mag.set_log_level(100)
    set_dtype(dtype_arg)

    mesh = mag.Mesh((nx, ny, nz), (500e-9/nx, 125e-9/ny, 3e-9))
    state = mag.State(mesh)
    state.material = {"Ms": 8e5, "A": 1.3e-11, "alpha": 0.02}

    demag = mag.DemagField()
    exchange = mag.ExchangeField()
    external = mag.ExternalField([-24.6e-3/mag.constants.mu_0,
                                  +4.3e-3/mag.constants.mu_0, 0.0])

    # s-state-like initial magnetization (skip the relax phase)
    state.m = state.Constant([0, 0, 0])
    state.m[1:-1, :, :, 0] = 1.0
    state.m[(-1, 0), :, :, 1] = 1.0
    mag.normalize(state.m)

    llg = mag.LLGSolver([demag, exchange, external])

    dt = 1e-11
    warmup = 20 if torch.cuda.is_available() else 5
    steps = 500 if torch.cuda.is_available() else 25

    for _ in range(warmup):
        llg.step(state, dt)

    reset_peak_memory()
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    t0 = time.perf_counter()
    for _ in range(steps):
        llg.step(state, dt)
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    elapsed = time.perf_counter() - t0

    print(json.dumps({
        "key": "sp4/%s/%dx%dx%d" % (dtype_arg, nx, ny, nz),
        "steps_per_s": round(steps / elapsed, 2),
        "m_avg": state.avg(state.m).tolist(),
        **peak_memory(),
    }))


if __name__ == "__main__":
    if sys.argv[1:2] == ["--variant"]:
        dtype_arg = sys.argv[2]
        run_variant(dtype_arg, *map(int, sys.argv[3:6]))
    else:
        from common import run_variants, out_dir_from_args
        out = out_dir_from_args()
        variants = [(d,) + s for s in SIZES for d in ("f64", "f32")]
        run_variants(os.path.abspath(__file__), variants, out, "sp4")
