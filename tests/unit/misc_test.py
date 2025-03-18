import pytest
from magnumnp import *
from helpers import *

def test_log_cumsum(simple_state):
    simple_state.x = torch.tensor(1.)
    log_func = lambda state: state.x
    cumsum_func = LogCumSum(log_func)
    assert cumsum_func(simple_state).cpu() == pytest.approx(1)
    assert cumsum_func(simple_state).cpu() == pytest.approx(2)
    assert cumsum_func(simple_state).cpu() == pytest.approx(3)

def test_log_diff():
    log_func = lambda state: state.x
    diff_func = LogDt(log_func)

    simple_state.t = 1
    simple_state.x = torch.tensor(1.)
    assert diff_func(simple_state).cpu() == pytest.approx(1)
    simple_state.t = 2
    simple_state.x = torch.tensor(1.)
    assert diff_func(simple_state).cpu() == pytest.approx(0)
    simple_state.t = 3
    simple_state.x = torch.tensor(2.)
    assert diff_func(simple_state).cpu() == pytest.approx(1)

def test_log_moving_average():
    log_func = lambda state: state.x
    avg_func = LogMovingAverage(log_func, N=3)

    simple_state.x = torch.tensor(1.)
    assert avg_func(simple_state) == pytest.approx(1.0)
    simple_state.x = torch.tensor(2.)
    assert avg_func(simple_state) == pytest.approx(1.5)
    simple_state.x = torch.tensor(3.)
    assert avg_func(simple_state) == pytest.approx(2.0)
