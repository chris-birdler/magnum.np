import pytest
from magnumnp import *

def test_log_cumsum():
    assert log_cumsum(x = 1) == pytest.approx(1)    
    assert log_cumsum(x = 1) == pytest.approx(2)    
    assert log_cumsum(x = 1) == pytest.approx(3)    

    assert log_cumsum(y = 1) == pytest.approx(1)    
    assert log_cumsum(y = 1) == pytest.approx(2)    
    assert log_cumsum(y = 1) == pytest.approx(3)    

def test_log_diff():
    assert log_diff(x = 1) == pytest.approx(1)    
    assert log_diff(x = 1) == pytest.approx(0)    
    assert log_diff(x = 2) == pytest.approx(1)    

    assert log_diff(y = 1) == pytest.approx(1)    
    assert log_diff(y = 1) == pytest.approx(0)    
    assert log_diff(y = 2) == pytest.approx(1)    



