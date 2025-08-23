import numpy as np
from cliffordlab.d_regular.hyperbolic import d_mul, d_inv, combine_idempotent, split_idempotent
from cliffordlab.d_regular.builders import from_FG

def test_inverse_roundtrip():
    a = (np.array([[2.0]]), np.array([[0.5]]))
    inva = d_inv(a)
    one = d_mul(a, inva)
    assert np.allclose(one[0], 1.0, atol=1e-6)
    assert np.allclose(one[1], 0.0, atol=1e-6)

def test_idempotent_roundtrip():
    p = np.array([[1.0, 2.0]])
    q = np.array([[0.3, -0.4]])
    F, G = split_idempotent(p, q)
    pp, qq = combine_idempotent(F, G)
    assert np.allclose(p, pp)
    assert np.allclose(q, qq)

def test_builder_FG():
    x = np.linspace(-1, 1, 11); y = np.linspace(-1, 1, 11)
    X, Y = np.meshgrid(x, y)
    F = lambda u: u**2
    G = lambda v: 2*v
    p, q = from_FG(F, G, X, Y)
    # spot check a point
    u = X[5,5] + Y[5,5]; v = X[5,5] - Y[5,5]
    p0 = 0.5 * (u**2 + 2*v)
    q0 = 0.5 * (u**2 - 2*v)
    assert np.allclose(p[5,5], p0) and np.allclose(q[5,5], q0)
