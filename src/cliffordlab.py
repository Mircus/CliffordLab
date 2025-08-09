from __future__ import annotations
from typing import Dict, Tuple, Iterable, Union
import math

Number = Union[float, int]
try:
    import sympy as sp
    Sym = sp.Expr
except Exception:
    sp = None
    Sym = Number

class Cl:
    def __init__(self, p: int, q: int):
        assert p >= 0 and q >= 0
        self.p, self.q, self.n = p, q, p+q
        self.sig = [1]*p + [-1]*q

    def basis_vector(self, i: int):
        assert 1 <= i <= self.n
        return MV.from_mask(self, 1, 1 << (i-1))

    def scalar(self, a: Union[Number, Sym]): return MV.from_mask(self, a, 0)
    def zero(self): return self.scalar(0)
    def vector(self, comps: Iterable[Union[Number, Sym]]):
        comps = list(comps); assert len(comps)==self.n
        out = self.zero()
        for i,c in enumerate(comps, start=1):
            if c: out = out + self.basis_vector(i)*c
        return out

    def gp_blades(self, a_mask: int, b_mask: int) -> Tuple[int, int]:
        sign = 1; k = 0
        for i in range(self.n):
            if (a_mask >> i) & 1:
                lower = b_mask & ((1 << i) - 1)
                k += bin(lower).count("1")
        if k & 1: sign = -sign
        c_mask = a_mask ^ b_mask
        overlap = a_mask & b_mask; i = 0
        while overlap:
            if overlap & 1: sign *= self.sig[i]
            overlap >>= 1; i += 1
        return sign, c_mask

class MV:
    def __init__(self, alg: Cl, data: Dict[int, Union[Number, Sym]] = None):
        self.A = alg; self.d = dict(data) if data else {}

    @staticmethod
    def from_mask(alg: Cl, coeff: Union[Number, Sym], mask: int):
        return MV(alg, {mask: coeff} if coeff else {})

    def sc(self): return self.d.get(0, 0)
    def grade_mask(self, m: int) -> int: return bin(m).count("1")
    def grade(self, k: int): return MV(self.A, {m:c for m,c in self.d.items() if self.grade_mask(m)==k})

    def __add__(self, other):
        self._check(other); out = dict(self.d)
        for m,c in other.d.items():
            out[m] = out.get(m,0) + c
            if out[m] == 0: del out[m]
        return MV(self.A, out)

    def __sub__(self, other):
        self._check(other); out = dict(self.d)
        for m,c in other.d.items():
            out[m] = out.get(m,0) - c
            if out[m] == 0: del out[m]
        return MV(self.A, out)

    def __mul__(self, other):
        if isinstance(other, MV):
            self._check(other); out = {}
            for am,ac in self.d.items():
                for bm,bc in other.d.items():
                    sgn, cm = self.A.gp_blades(am, bm)
                    val = ac*bc*sgn
                    out[cm] = out.get(cm,0) + val
                    if out[cm]==0: del out[cm]
            return MV(self.A, out)
        else:
            if other == 0: return self.A.zero()
            return MV(self.A, {m:c*other for m,c in self.d.items()})

    __rmul__ = __mul__
    def __truediv__(self, s): return MV(self.A, {m:c/s for m,c in self.d.items()})
    def __neg__(self): return MV(self.A, {m:-c for m,c in self.d.items()})
    def _check(self, other): assert self.A is other.A, "Algebras differ"

    def reverse(self):
        out = {}
        for m,c in self.d.items():
            k = bin(m).count("1")
            sgn = -1 if ((k*(k-1)//2)%2) else 1
            out[m] = sgn*c
        return MV(self.A, out)

    def grade_involution(self):
        out = {}
        for m,c in self.d.items():
            k = bin(m).count("1")
            out[m] = (-1 if k%2 else 1)*c
        return MV(self.A, out)

    def conjugate(self): return self.grade_involution().reverse()

    def inner(self, other):
        self._check(other); out = self.A.zero().d
        for am,ac in self.d.items():
            a = bin(am).count("1")
            for bm,bc in other.d.items():
                sgn, cm = self.A.gp_blades(am, bm)
                k = bin(cm).count("1"); b = bin(bm).count("1")
                if k == abs(a-b):
                    out[cm] = out.get(cm,0) + ac*bc*sgn
        return MV(self.A, out)

    def wedge(self, other):
        self._check(other); out = self.A.zero().d
        for am,ac in self.d.items():
            a = bin(am).count("1")
            for bm,bc in other.d.items():
                sgn, cm = self.A.gp_blades(am, bm)
                k = bin(cm).count("1"); b = bin(bm).count("1")
                if k == a + b:
                    out[cm] = out.get(cm,0) + ac*bc*sgn
        return MV(self.A, out)

    def __repr__(self):
        if not self.d: return "0"
        terms=[]
        for m in sorted(self.d.keys()):
            c=self.d[m]
            if m==0: terms.append(f"{c}")
            else:
                idxs=[str(i+1) for i in range(self.A.n) if (m>>i)&1]
                terms.append(f"{c}*e{''.join(idxs)}")
        return " + ".join(terms)

def unit_J_from_axis(A: Cl, axis: int, sign: int) -> MV:
    assert 1 <= axis <= A.n and sign in (-1,+1)
    e = A.basis_vector(axis); sq = A.sig[axis-1]
    if sq != sign: raise ValueError(f"e{axis}^2={sq} but requested J^2={sign}")
    return e

def slice_mul(a, b, sign: int):
    x1,y1=a; x2,y2=b
    return (x1*x2 + sign*y1*y2, x1*y2 + x2*y1)

def slice_inv(a, sign: int):
    x,y = a; den = x*x - sign*y*y
    if den == 0: raise ZeroDivisionError("Non-invertible on slice")
    return (x/den, -y/den)

def slice_exp(a, sign: int, N: int=48):
    out=(1,0); term=(1,0)
    for n in range(1,N+1):
        term = slice_mul(term, a, sign)
        out = (out[0] + term[0]/math.factorial(n), out[1] + term[1]/math.factorial(n))
    return out

def CR_residuals(fAB, x: float, y: float, sign: int, h: float=1e-5):
    A=lambda X,Y: fAB(X,Y,sign)[0]; B=lambda X,Y: fAB(X,Y,sign)[1]
    dAx=(A(x+h,y)-A(x-h,y))/(2*h); dAy=(A(x,y+h)-A(x,y-h))/(2*h)
    dBx=(B(x+h,y)-B(x-h,y))/(2*h); dBy=(B(x,y+h)-B(x,y-h))/(2*h)
    return (dAx - dBy, dAy + dBx) if sign==-1 else (dAx - dBy, dAy - dBx)

def cauchy_reconstruct(fAB, z0, sign: int, radius: float=0.35, N: int=800):
    x0,y0=z0; acc=(0.0,0.0)
    import numpy as np
    if sign==-1:
        ts=np.linspace(0,2*math.pi,N,endpoint=False)
        for k in range(N):
            t1=ts[k]; t2=ts[(k+1)%N]; tm=0.5*(t1+t2)
            z=(x0+radius*math.cos(tm), y0+radius*math.sin(tm))
            ker=slice_inv((z[0]-x0,z[1]-y0),sign)
            val=fAB(z[0],z[1],sign); a,b = slice_mul(ker,val,sign)
            dt=t2-t1; ds=radius*dt
            acc=(acc[0]+(-b*sign*ds), acc[1]+(-a*ds))
        return (acc[0]/(2*math.pi), acc[1]/(2*math.pi))
    else:
        r=radius; steps=N//4
        pts=[((x0-r,y0+r),(x0+r,y0+r)),((x0+r,y0+r),(x0+r,y0-r)),
             ((x0+r,y0-r),(x0-r,y0-r)),((x0-r,y0-r),(x0-r,y0+r))]
        for (xa,ya),(xb,yb) in pts:
            xs=np.linspace(xa,xb,steps,endpoint=False); ys=np.linspace(ya,yb,steps,endpoint=False)
            for i in range(steps):
                xm=0.5*(xs[i]+xs[(i+1)%steps]); ym=0.5*(ys[i]+ys[(i+1)%steps])
                z=(xm,ym); ker=slice_inv((z[0]-x0,z[1]-y0),sign)
                val=fAB(z[0],z[1],sign); a,b = slice_mul(ker,val,sign)
                dx=(xb-xa)/steps; dy=(yb-ya)/steps; ds=math.hypot(dx,dy)
                acc=(acc[0]+(b*sign*ds), acc[1]+(a*ds))
        return (acc[0]/(2*math*pi), acc[1]/(2*math*pi))

def demo_basic():
    A=Cl(3,1); e4=A.basis_vector(4); e1=A.basis_vector(1)
    def fAB(x,y,sign): return slice_exp((x,y),sign,N=40)
    r_m=CR_residuals(fAB,0.2,0.3,-1); r_p=CR_residuals(fAB,0.2,0.3,+1)
    rec_m=cauchy_reconstruct(fAB,(0.2,0.3),-1); rec_p=cauchy_reconstruct(fAB,(0.2,0.3),+1)
    return {"CR_residuals":{"elliptic":r_m,"hyperbolic":r_p},
            "Cauchy_reconstruction":{"elliptic":{"recon":rec_m},"hyperbolic":{"recon":rec_p}}}
