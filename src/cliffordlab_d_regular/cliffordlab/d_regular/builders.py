import numpy as np
from .hyperbolic import combine_idempotent, u_of, v_of

__all__ = ["from_FG", "poly_FG", "exp_linear", "sqrt_abs_regularized"]

def from_FG(F, G, x, y):
    """
    Build a D-holomorphic function f = F(u) e_+ + G(v) e_-, sampled on grid (x,y).
    F, G : callables on R -> R
    Returns (p,q) arrays with f = p + j q.
    """
    u = u_of(x, y); v = v_of(x, y)
    Fu = F(u); Gv = G(v)
    return combine_idempotent(Fu, Gv)

def poly_FG(coeffs_F, coeffs_G, x, y):
    """
    Polynomial F(u) = sum c_k u^k, G(v) = sum d_k v^k.
    coeffs_*: iterable [c0, c1, ..., cn].
    """
    u = u_of(x, y); v = v_of(x, y)
    Fu = sum(c * u**k for k, c in enumerate(coeffs_F))
    Gv = sum(d * v**k for k, d in enumerate(coeffs_G))
    return combine_idempotent(Fu, Gv)

def exp_linear(a_plus, a_minus, x, y):
    """
    F(u)=exp(a_plus * u), G(v)=exp(a_minus * v).
    """
    u = u_of(x, y); v = v_of(x, y)
    Fu = np.exp(a_plus * u); Gv = np.exp(a_minus * v)
    return combine_idempotent(Fu, Gv)

def sqrt_abs_regularized(eps=1e-6):
    """
    Return a function H(t) ~ sign(t) * sqrt(|t|) with small smoothing near 0.
    Useful to create pieces with controlled derivative jump along characteristics.
    """
    def H(t):
        t = np.asarray(t)
        return np.sign(t) * np.sqrt(np.sqrt(t**2 + eps**2))
    return H
