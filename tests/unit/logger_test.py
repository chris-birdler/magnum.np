import pytest
import torch
from magnumnp import *
from helpers import *
import numpy as np

def test_log(simple_state, tmpdir):
    p = str(tmpdir)
    logger = Logger(p, scalars = ["t", "m"], fields = ["m"])
    logger << simple_state

def test_is_resumable(simple_state, tmpdir):
    p = str(tmpdir)
    logger = Logger(p, scalars = ["t", "m"], fields = ["m"])
    assert logger.is_resumable() == False
    assert logger.loggers["fields"].last_recorded_step() == None
    logger << simple_state
    assert logger.is_resumable() == True
    assert logger.loggers["fields"].last_recorded_step() == 0
    logger << simple_state
    assert logger.is_resumable() == True
    assert logger.loggers["fields"].last_recorded_step() == 1

def test_resume(simple_state, tmpdir):
    p = str(tmpdir / "log")

    logger = Logger(p, scalars = ["t", "m"], fields = ["m"], fields_every = 8)
    for t in torch.arange(0., 2.8, 0.1):
        simple_state.t = float(t)
        simple_state.m = simple_state.Constant([t, 0, 0])
        logger << simple_state

    del logger

    rlogger = Logger(p, scalars = ["t", "m"], fields = ["m"], fields_every = 8)
    assert rlogger.loggers["scalars"].resumable_step() == 28
    assert rlogger.loggers["fields"].resumable_step() == 32
    rlogger.resume(simple_state)

    torch.testing.assert_close(simple_state.avg(simple_state.m), torch.tensor([2.4,0.,0.]))
    assert simple_state.t == pytest.approx(2.4)
    for t in torch.arange(0., 0.8, 0.1):
        simple_state.t = float(t)
        simple_state.m = simple_state.Constant([t, 0, 0])
        rlogger << simple_state

def test_nonequi(tmpdir):
    n  = (20, 4, 4)
    dx = (5., 3., np.arange(n[2]) + 1.)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.m = state.Constant([1,0,0])

    p = str(tmpdir)
    logger = Logger(p, scalars = ["t", "m"], fields = ["m"])
    logger << state
