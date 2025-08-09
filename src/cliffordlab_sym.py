from __future__ import annotations
from typing import Tuple, Union
import sympy as sp

Number = Union[float, int]
Sym = sp.Expr

def slice_mul(a: Tuple[Sym, Sym], b: Tuple[Sym, Sym], sign: int) -> Tuple[Sym, Sym]:
    x1,y1=a; x2,y2=b
    return (x1*x2 + sign*y1*y2, x1*y2 + x2*y1)

def slice_exp(a: Tuple[Sym, Sym], sign: int, N: int=60) -> Tuple[Sym, Sym]:
    out=(sp.Integer(1), sp.Integer(0)); term=(sp.Integer(1), sp.Integer(0))
    for n in range(1,N+1):
        term = slice_mul(term, a, sign)
        out = (out[0] + term[0]/sp.factorial(n), out[1] + term[1]/sp.factorial(n))
    return sp.simplify(out[0]), sp.simplify(out[1])

def CR_residuals_symbolic(A: Sym, B: Sym, x: Sym, y: Sym, sign: int) -> Tuple[Sym, Sym]:
    dAx, dAy = sp.diff(A,x), sp.diff(A,y)
    dBx, dBy = sp.diff(B,x), sp.diff(B,y)
    if sign==-1: return (sp.simplify(dAx - dBy), sp.simplify(dAy + dBx))
    else:        return (sp.simplify(dAx - dBy), sp.simplify(dAy - dBx))

def example_exp_symbolic(sign: int, N: int=40):
    x,y = sp.symbols('x y', real=True)
    A,B = slice_exp((x,y), sign, N=N)
    return (x,y), (A,B), CR_residuals_symbolic(A,B,x,y,sign)

def example_exp_closed_forms():
    x,y = sp.symbols('x y', real=True)
    Ae,Be = sp.exp(x)*sp.cos(y), sp.exp(x)*sp.sin(y)
    Ah,Bh = sp.exp(x)*sp.cosh(y), sp.exp(x)*sp.sinh(y)
    r1e,r2e = CR_residuals_symbolic(Ae,Be,x,y,-1)
    r1h,r2h = CR_residuals_symbolic(Ah,Bh,x,y,+1)
    return (x,y),(Ae,Be,r1e,r2e),(Ah,Bh,r1h,r2h)

if __name__ == "__main__":
    (x,y),(A,B),(r1,r2) = example_exp_symbolic(-1,60)
    print("Elliptic CR residuals:", sp.simplify(r1), sp.simplify(r2))
    (x,y),(A,B),(r1,r2) = example_exp_symbolic(+1,60)
    print("Hyperbolic CR residuals:", sp.simplify(r1), sp.simplify(r2))
