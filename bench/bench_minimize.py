"""MinimizerBB iteration throughput (f32/f64)."""
import json
import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def run_variant(dtype_arg):
    import torch
    import magnumnp as mag
    from common import reset_peak_memory, peak_memory, set_dtype

    mag.set_log_level(100)
    set_dtype(dtype_arg)

    mesh = mag.Mesh((200, 50, 1), (2.5e-9, 2.5e-9, 3e-9))
    state = mag.State(mesh)
    state.material = {"Ms": 8e5, "A": 1.3e-11, "alpha": 0.02}

    demag = mag.DemagField()
    exchange = mag.ExchangeField()

    state.m = state.Constant([0, 0, 0])
    state.m[1:-1, :, :, 0] = 1.0
    state.m[(-1, 0), :, :, 1] = 1.0
    mag.normalize(state.m)

    minimizer = mag.MinimizerBB([demag, exchange])
    minimizer.minimize(state, maxiter=20) # warmup (includes kernel init + compile)

    iters = 200 if torch.cuda.is_available() else 50
    reset_peak_memory()
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    t0 = time.perf_counter()
    minimizer.minimize(state, maxiter=iters, dm_tol=0.) # dm_tol=0 -> never converges early
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    elapsed = time.perf_counter() - t0

    print(json.dumps({
        "key": "minimize/%s" % dtype_arg,
        "iters_per_s": round(iters / elapsed, 2),
        **peak_memory(),
    }))


if __name__ == "__main__":
    if sys.argv[1:2] == ["--variant"]:
        run_variant(sys.argv[2])
    else:
        from common import run_variants, out_dir_from_args
        out = out_dir_from_args()
        run_variants(os.path.abspath(__file__), [("f64",), ("f32",)], out, "minimize")
