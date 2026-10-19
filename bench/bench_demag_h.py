"""Benchmark DemagField.h: per-component FFT loop vs batched single-FFT, f32/f64.

Decides the batched-FFT question: adopt only if it wins >=10% at BOTH precisions.
"""
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SIZES = [(128, 128, 8), (256, 256, 16)]


def run_variant(dtype_arg, impl, nx, ny, nz):
    import torch
    import magnumnp as mag
    from common import time_fn, reset_peak_memory, peak_memory, set_dtype

    mag.set_log_level(100)
    set_dtype(dtype_arg)

    def h_batched(self, state):
        if not hasattr(self, "_N"):
            self._N = self._init_N(state)
        dim = [i for i in range(3) if state.mesh.n[i] > 1]
        shape = self._shape(state)
        s = [shape[i] for i in dim]

        M = state.material["Ms"] * state.m
        m_fft = torch.fft.rfftn(M, dim=dim, s=s)
        mx, my, mz = m_fft.unbind(-1)
        N = self._N
        hf = torch.stack([N[0][0]*mx + N[0][1]*my + N[0][2]*mz,
                          N[1][0]*mx + N[1][1]*my + N[1][2]*mz,
                          N[2][0]*mx + N[2][1]*my + N[2][2]*mz], dim=-1)
        del m_fft, mx, my, mz
        h = torch.fft.irfftn(hf, dim=dim)
        return h[:state.mesh.n[0], :state.mesh.n[1], :state.mesh.n[2]].contiguous()

    if impl == "batched":
        mag.DemagField.h = h_batched

    mesh = mag.Mesh((nx, ny, nz), (5e-9, 5e-9, 5e-9))
    state = mag.State(mesh)
    state.material = {"Ms": 8e5}
    state.m = state.RandM()

    demag = mag.DemagField()
    demag.h(state) # includes kernel init

    reset_peak_memory()
    h_ms = time_fn(lambda: demag.h(state), warmup=5, iters=50)
    checksum = float(demag.h(state).double().abs().sum())

    print(json.dumps({
        "key": "demag_h/%s/%s/%dx%dx%d" % (dtype_arg, impl, nx, ny, nz),
        "h_ms": round(h_ms, 3),
        "checksum": checksum,
        **peak_memory(),
    }))


if __name__ == "__main__":
    if sys.argv[1:2] == ["--variant"]:
        dtype_arg, impl, nx, ny, nz = sys.argv[2:7]
        run_variant(dtype_arg, impl, int(nx), int(ny), int(nz))
    else:
        from common import run_variants, out_dir_from_args
        out = out_dir_from_args()
        variants = [(d, i) + s for s in SIZES for d in ("f64", "f32") for i in ("loop", "batched")]
        run_variants(os.path.abspath(__file__), variants, out, "demag_h")
