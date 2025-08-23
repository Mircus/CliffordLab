import numpy as np

__all__ = [
    "u_of", "v_of", "split_idempotent", "combine_idempotent",
    "d_add", "d_sub", "d_mul", "d_inv", "d_conj",
    "d_exp", "d_log", "d_sin", "d_cos", "d_sinh", "d_cosh"
]

def u_of(x, y):
    """u = x + y"""
    return np.asarray(x) + np.asarray(y)

def v_of(x, y):
    """v = x - y"""
    return np.asarray(x) - np.asarray(y)

def split_idempotent(p, q):
    """
    Split f = p + j q into idempotent components:
    f = F(u) e_+ + G(v) e_- where
    F = p + q, G = p - q  (as functions sampled at (x,y))
    Returns (F_component, G_component). Note: these are evaluations, not functions.
    """
    p = np.asarray(p); q = np.asarray(q)
    F = p + q
    G = p - q
    return F, G

def combine_idempotent(F, G):
    """
    Combine idempotent components back to (p, q):
    p = (F + G)/2, q = (F - G)/2
    """
    F = np.asarray(F); G = np.asarray(G)
    p = 0.5 * (F + G)
    q = 0.5 * (F - G)
    return p, q

# Basic D-arithmetic for scalars/arrays: f = p + j q
def d_add(a, b):
    ap, aq = a
    bp, bq = b
    return (np.asarray(ap) + np.asarray(bp), np.asarray(aq) + np.asarray(bq))

def d_sub(a, b):
    ap, aq = a
    bp, bq = b
    return (np.asarray(ap) - np.asarray(bp), np.asarray(aq) - np.asarray(bq))

def d_mul(a, b):
    """
    (p + j q)(r + j s) = (pr + qs) + j(ps + qr) since j^2 = +1
    """
    p, q = a; r, s = b
    p = np.asarray(p); q = np.asarray(q); r = np.asarray(r); s = np.asarray(s)
    return (p*r + q*s, p*s + q*r)

def d_conj(a):
    """Conjugation j -> -j: (p, q) -> (p, -q)."""
    p, q = a
    return (np.asarray(p), -np.asarray(q))

def d_inv(a):
    """
    Inverse when not on zero divisors: (p + j q)^(-1) = (p - j q) / (p^2 - q^2).
    Returns (pinv, qinv). If denominator is zero, returns (inf, inf) by numpy convention.
    """
    p, q = a
    denom = np.asarray(p)**2 - np.asarray(q)**2
    return (np.asarray(p)/denom, -np.asarray(q)/denom)

# Special functions via standard identities
def d_exp(a):
    """exp(p + j q) = exp(p)*(cosh(q) + j*sinh(q))"""
    p, q = a
    ep = np.exp(np.asarray(p))
    return (ep * np.cosh(q), ep * np.sinh(q))

def d_log(a):
    """
    log(p + j q) in a region avoiding zero divisors.
    Use idempotent chart: f = F e_+ + G e_- with F = p+q, G = p-q
    Then log f = (log F) e_+ + (log G) e_- => back to (p,q).
    """
    F, G = split_idempotent(*a)
    # Guard: avoid <=0 for real log; users should choose branches as needed.
    logF = np.log(F)
    logG = np.log(G)
    return combine_idempotent(logF, logG)

def d_sinh(a):
    """sinh(p + j q) = sinh p * cosh q + j (cosh p * sinh q)"""
    p, q = a
    return (np.sinh(p)*np.cosh(q), np.cosh(p)*np.sinh(q))

def d_cosh(a):
    """cosh(p + j q) = cosh p * cosh q + j (sinh p * sinh q)"""
    p, q = a
    return (np.cosh(p)*np.cosh(q), np.sinh(p)*np.sinh(q))

def d_sin(a):
    """sin(p + j q) = sin p * cosh q + j (cos p * sinh q)"""
    p, q = a
    return (np.sin(p)*np.cosh(q), np.cos(p)*np.sinh(q))

def d_cos(a):
    """cos(p + j q) = cos p * cosh q - j (sin p * sinh q)"""
    p, q = a
    return (np.cos(p)*np.cosh(q), -np.sin(p)*np.sinh(q))
