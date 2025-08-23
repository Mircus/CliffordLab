import numpy as np
from .hyperbolic import u_of, v_of, combine_idempotent

__all__ = ["pv_u", "pv_v", "kernel_eval"]

def pv_u(x, y, alpha=0.0, eps=1e-6):
    """
    Approximate principal value 1/(u - alpha) by symmetric epsilon-exclusion:
    1/(u-a) ~ (u-a) / ((u-a)^2 + eps^2)
    """
    u = u_of(x, y)
    return (u - alpha) / ((u - alpha)**2 + eps**2)

def pv_v(x, y, beta=0.0, eps=1e-6):
    v = v_of(x, y)
    return (v - beta) / ((v - beta)**2 + eps**2)

def kernel_eval(x, y, alpha=0.0, beta=0.0, eps=1e-6):
    """
    Idempotent-decomposed hyperfunction-inspired kernel:
    K = (1/(u-alpha)) e_+ + (1/(v-beta)) e_-  (regularized principal value).
    Returns (p, q).
    """
    Fu = pv_u(x, y, alpha=alpha, eps=eps)
    Gv = pv_v(x, y, beta=beta, eps=eps)
    return combine_idempotent(Fu, Gv)
