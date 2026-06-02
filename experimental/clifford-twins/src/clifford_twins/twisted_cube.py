from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


def bit(a: int, i: int) -> int:
    return (a >> i) & 1


def popcount(a: int) -> int:
    return a.bit_count()


def dot_mod2(a: int, b: int) -> int:
    return popcount(a & b) % 2


def vectors(n: int) -> list[int]:
    return list(range(1 << n))


SignRule = Callable[[int, int], int]


@dataclass(frozen=True)
class Twist:
    """Boolean-cube sign-twist algebra.

    Basis indexed by a in (Z2)^n:

        u_a u_b = F(a,b) u_{a xor b}.

    The associator factor is:

        Phi(a,b,c)=F(a,b)F(a+b,c)/(F(b,c)F(a,b+c)).

    This supports the article's "twisted hypercube" viewpoint, but the main
    repo focus is the invariant Cl(3,0)=H⊕HI active-volume construction.
    """

    n: int
    F_rule: SignRule

    def F(self, a: int, b: int) -> int:
        value = self.F_rule(a, b)
        if value not in (-1, 1):
            raise ValueError("A sign twist must return ±1.")
        return value

    def multiply_basis(self, a: int, b: int) -> tuple[int, int]:
        return self.F(a, b), a ^ b

    def associator_factor(self, a: int, b: int, c: int) -> int:
        numerator = self.F(a, b) * self.F(a ^ b, c)
        denominator = self.F(b, c) * self.F(a, b ^ c)
        return numerator // denominator

    def is_associative(self) -> bool:
        for a in vectors(self.n):
            for b in vectors(self.n):
                for c in vectors(self.n):
                    if self.associator_factor(a, b, c) != 1:
                        return False
        return True
