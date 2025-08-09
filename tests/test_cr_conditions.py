import math
from cliffordlab import CR_residuals, slice_exp

def test_exp_cr_elliptic():
    def fAB(x,y,s): return slice_exp((x,y), s, N=30)
    r1,r2 = CR_residuals(fAB, 0.2, 0.3, -1)
    assert abs(r1) < 1e-3 and abs(r2) < 1e-3

def test_exp_cr_hyperbolic():
    def fAB(x,y,s): return slice_exp((x,y), s, N=30)
    r1,r2 = CR_residuals(fAB, 0.2, 0.3, +1)
    assert abs(r1) < 1e-3 and abs(r2) < 1e-3
