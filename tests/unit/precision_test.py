import json
import os
import subprocess
import sys
import textwrap

import pytest

# set_precision changes the global default dtype and the pytest process has
# already created float64 meshes, so each precision variant runs in its own
# subprocess and reports its results as json.
SCRIPT = textwrap.dedent("""
    import json, sys
    from magnumnp import *
    import torch

    set_log_level(100)
    set_precision(sys.argv[1])

    n  = (25, 10, 1)
    mesh = Mesh(n, (5e-9, 5e-9, 3e-9))
    state = State(mesh)
    state.material = {"Ms": 8e5, "A": 1.3e-11, "alpha": 1.0}

    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField([-24.6e-3/constants.mu_0, +4.3e-3/constants.mu_0, 0.0])

    # s-state-like initial magnetization
    state.m = state.Constant([0,0,0])
    state.m[1:-1,:,:,0]   = 1.0
    state.m[(-1,0),:,:,1] = 1.0

    llg = LLGSolver([demag, exchange, external])
    for i in range(5):
        llg.step(state, 1e-11)

    m_norm = torch.linalg.norm(state.m, dim=-1)
    print(json.dumps({
        "dtype": str(state.dtype),
        "m_avg": state.avg(state.m).tolist(),
        "m_norm_min": float(m_norm.min()),
        "m_norm_max": float(m_norm.max()),
    }))
""")


def run_variant(precision):
    result = subprocess.run([sys.executable, "-c", SCRIPT, precision],
                            capture_output=True, text=True, timeout=600,
                            env={**os.environ, "CUDA_DEVICE": os.environ.get("CUDA_DEVICE", "-1")})
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout.strip().splitlines()[-1])


def test_single_precision():
    single = run_variant("single")
    double = run_variant("double")

    assert single["dtype"] == "torch.float32"
    assert double["dtype"] == "torch.float64"

    # normalization survives single precision
    assert single["m_norm_min"] == pytest.approx(1., abs=1e-5)
    assert single["m_norm_max"] == pytest.approx(1., abs=1e-5)

    # trajectory agrees with the double precision run
    for s, d in zip(single["m_avg"], double["m_avg"]):
        assert s == pytest.approx(d, abs=5e-3)


def test_set_precision_after_mesh_raises():
    script = textwrap.dedent("""
        from magnumnp import *
        set_log_level(100)
        mesh = Mesh((4, 4, 1), (1e-9, 1e-9, 1e-9))
        try:
            set_precision("single")
        except RuntimeError:
            print("RAISED")
        set_precision("single", force=True) # escape hatch still works
        print("FORCED")
    """)
    result = subprocess.run([sys.executable, "-c", script],
                            capture_output=True, text=True, timeout=300,
                            env={**os.environ, "CUDA_DEVICE": os.environ.get("CUDA_DEVICE", "-1")})
    assert result.returncode == 0, result.stderr
    assert "RAISED" in result.stdout
    assert "FORCED" in result.stdout
