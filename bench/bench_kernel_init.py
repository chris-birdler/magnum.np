"""Benchmark demag kernel initialization: wall time + peak (GPU) memory.

Headline metric for the chunked CPU-resident kernel setup: the device peak
during init should be close to the final kernel size instead of 10-20x.
"""
import json
import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SIZES = [(128, 128, 8), (256, 256, 16), (512, 512, 4)]


def run_variant(nx, ny, nz):
    import magnumnp as mag
    from common import reset_peak_memory, peak_memory

    mag.set_log_level(100)

    mesh = mag.Mesh((nx, ny, nz), (5e-9, 5e-9, 5e-9))
    state = mag.State(mesh)
    state.material = {"Ms": 8e5}

    reset_peak_memory()
    t0 = time.perf_counter()
    N = mag.DemagField()._init_N(state)
    init_s = time.perf_counter() - t0

    kernel_mb = sum(c.numel() * c.element_size() for row in N for c in row) * 6 // 9 / 1e6

    print(json.dumps({
        "key": "kernel_init/%dx%dx%d" % (nx, ny, nz),
        "init_s": round(init_s, 2),
        "kernel_mb": round(kernel_mb, 1),
        **peak_memory(),
    }))


if __name__ == "__main__":
    if sys.argv[1:2] == ["--variant"]:
        run_variant(*map(int, sys.argv[2:5]))
    else:
        from common import run_variants, out_dir_from_args
        out = out_dir_from_args()
        run_variants(os.path.abspath(__file__), SIZES, out, "kernel_init")
