from __future__ import annotations

from dataclasses import dataclass
from math import isclose, sqrt


@dataclass(frozen=True)
class Quaternion:
    """A small dependency-free quaternion class.

    Coordinates are ordered as:

        r + x*i + y*j + z*k.

    In the article's language this models:

        H = Cl^+(3,0) ≅ H,

    the quaternionic rotor/even Clifford core.
    """

    r: float = 0.0
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    def __add__(self, other: Quaternion) -> Quaternion:
        return Quaternion(
            self.r + other.r,
            self.x + other.x,
            self.y + other.y,
            self.z + other.z,
        )

    def __sub__(self, other: Quaternion) -> Quaternion:
        return Quaternion(
            self.r - other.r,
            self.x - other.x,
            self.y - other.y,
            self.z - other.z,
        )

    def __neg__(self) -> Quaternion:
        return Quaternion(-self.r, -self.x, -self.y, -self.z)

    def __mul__(self, other: Quaternion | float | int) -> Quaternion:
        if isinstance(other, (float, int)):
            return Quaternion(
                self.r * other,
                self.x * other,
                self.y * other,
                self.z * other,
            )

        a, b, c, d = self.r, self.x, self.y, self.z
        e, f, g, h = other.r, other.x, other.y, other.z

        return Quaternion(
            a * e - b * f - c * g - d * h,
            a * f + b * e + c * h - d * g,
            a * g - b * h + c * e + d * f,
            a * h + b * g - c * f + d * e,
        )

    def __rmul__(self, scalar: float | int) -> Quaternion:
        return self * scalar

    def conj(self) -> Quaternion:
        """Quaternionic conjugation.

        In Cl^+(3,0), this is Clifford reversion restricted to the even core.
        """
        return Quaternion(self.r, -self.x, -self.y, -self.z)

    def norm2(self) -> float:
        return self.r * self.r + self.x * self.x + self.y * self.y + self.z * self.z

    def norm(self) -> float:
        return sqrt(self.norm2())

    def approx_eq(self, other: Quaternion, tol: float = 1e-9) -> bool:
        return (
            isclose(self.r, other.r, abs_tol=tol)
            and isclose(self.x, other.x, abs_tol=tol)
            and isclose(self.y, other.y, abs_tol=tol)
            and isclose(self.z, other.z, abs_tol=tol)
        )

    def as_tuple(self) -> tuple[float, float, float, float]:
        return self.r, self.x, self.y, self.z

    def __repr__(self) -> str:
        return f"Quaternion({self.r:g}, {self.x:g}, {self.y:g}, {self.z:g})"
