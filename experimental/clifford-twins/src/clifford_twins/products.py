from __future__ import annotations

from collections.abc import Callable

from .cl3 import Cl3Element

Product = Callable[[Cl3Element, Cl3Element], Cl3Element]


def clifford_product(x: Cl3Element, y: Cl3Element) -> Cl3Element:
    """Passive-volume product: the Cl(3,0) product in H ⊕ H I.

    I is central and I^2 = -1:

        (a+bI)(c+dI) = (ac - bd) + (ad + bc)I.

    This is associative.
    """
    a, b = x.a, x.b
    c, d = y.a, y.b
    return Cl3Element(a * c - b * d, a * d + b * c)


def octonion_twin_product(x: Cl3Element, y: Cl3Element) -> Cl3Element:
    """Active-volume product: the octonion twin product on H ⊕ H ell.

    The active volume element obeys:

        ell a = conjugate(a) ell.

    We use the Cayley-Dickson convention:

        (a+b ell)(c+d ell)
        =
        (ac - conjugate(d)b) + (da + b conjugate(c))ell.

    This is the octonion product written on the same 8D vector space.
    """
    a, b = x.a, x.b
    c, d = y.a, y.b
    return Cl3Element(a * c - d.conj() * b, d * a + b * c.conj())


def associator(product: Product, x: Cl3Element, y: Cl3Element, z: Cl3Element) -> Cl3Element:
    """Associator [x,y,z] = (xy)z - x(yz)."""
    return product(product(x, y), z) - product(x, product(y, z))


def commutator(product: Product, x: Cl3Element, y: Cl3Element) -> Cl3Element:
    """Commutator [x,y] = xy - yx."""
    return product(x, y) - product(y, x)
