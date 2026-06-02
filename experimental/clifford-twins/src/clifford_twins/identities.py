from __future__ import annotations

from itertools import product as cartesian_product

from .basis import basis
from .cl3 import Cl3Element
from .products import Product, associator


def is_associative_on_basis(product: Product, tol: float = 1e-9) -> bool:
    """Check associativity on basis triples."""
    for x, y, z in cartesian_product(basis(), repeat=3):
        if not associator(product, x, y, z).is_zero(tol=tol):
            return False
    return True


def is_alternative_on_basis(product: Product, tol: float = 1e-9) -> bool:
    """Check alternativity on basis pairs:

        [x,x,y]=0 and [y,x,x]=0.
    """
    for x, y in cartesian_product(basis(), repeat=2):
        if not associator(product, x, x, y).is_zero(tol=tol):
            return False
        if not associator(product, y, x, x).is_zero(tol=tol):
            return False
    return True


def associator_table(product: Product) -> dict[tuple[int, int, int], Cl3Element]:
    """Return all nonzero basis associators.

    Indices refer to basis() ordering:
        0=1, 1=i, 2=j, 3=k, 4=ell, 5=iell, 6=jell, 7=kell.
    """
    B = basis()
    out: dict[tuple[int, int, int], Cl3Element] = {}
    for i, x in enumerate(B):
        for j, y in enumerate(B):
            for k, z in enumerate(B):
                a = associator(product, x, y, z)
                if not a.is_zero():
                    out[(i, j, k)] = a
    return out
