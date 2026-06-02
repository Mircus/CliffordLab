from __future__ import annotations

from .cl3 import Cl3Element
from .products import octonion_twin_product


def real_part(x: Cl3Element) -> float:
    """Real part of an octonion/Cl3Element."""
    return x.real()


def inner_product(x: Cl3Element, y: Cl3Element) -> float:
    """Euclidean inner product on the 8D vector space.

    This is the coordinate dot product on H⊕H.
    """
    return sum(a * b for a, b in zip(x.as_tuple(), y.as_tuple()))


def triality_form(x: Cl3Element, psi: Cl3Element, phi: Cl3Element) -> float:
    """Octonionic triality trilinear form.

    Schematically:

        t(x, psi, phi) = Re(x * (psi * phi)),

    using the octonion twin product.

    Note:
    Octonions are nonassociative, so this function intentionally fixes
    the bracketing x*(psi*phi). Other equivalent formulations may use
    related cyclic real-part identities.
    """
    return real_part(octonion_twin_product(x, octonion_twin_product(psi, phi)))
