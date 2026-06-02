from __future__ import annotations

from dataclasses import dataclass
from math import isclose

from .quaternion import Quaternion


@dataclass(frozen=True)
class Cl3Element:
    """Element of the shared 8D vector space H ⊕ H I.

    Write elements as:

        a + b I,

    where a,b in H = Cl^+(3,0) ≅ quaternions.

    This class is only the vector-space object. It deliberately does NOT
    choose a multiplication. Multiplications live in products.py.
    """

    a: Quaternion = Quaternion()
    b: Quaternion = Quaternion()

    def __add__(self, other: Cl3Element) -> Cl3Element:
        return Cl3Element(self.a + other.a, self.b + other.b)

    def __sub__(self, other: Cl3Element) -> Cl3Element:
        return Cl3Element(self.a - other.a, self.b - other.b)

    def __neg__(self) -> Cl3Element:
        return Cl3Element(-self.a, -self.b)

    def scale(self, scalar: float) -> Cl3Element:
        return Cl3Element(self.a * scalar, self.b * scalar)

    def conj_octonion(self) -> Cl3Element:
        """Cayley-Dickson conjugation on H ⊕ H ell."""
        return Cl3Element(self.a.conj(), -self.b)

    def real(self) -> float:
        return self.a.r

    def norm2(self) -> float:
        return self.a.norm2() + self.b.norm2()

    def approx_eq(self, other: Cl3Element, tol: float = 1e-9) -> bool:
        return self.a.approx_eq(other.a, tol=tol) and self.b.approx_eq(other.b, tol=tol)

    def is_zero(self, tol: float = 1e-9) -> bool:
        return self.approx_eq(Cl3Element(), tol=tol)

    def as_tuple(self) -> tuple[float, ...]:
        return (*self.a.as_tuple(), *self.b.as_tuple())

    def __repr__(self) -> str:
        return f"Cl3Element(a={self.a!r}, b={self.b!r})"
