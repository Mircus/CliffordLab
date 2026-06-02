"""Clifford Twins: active-volume calculus for Cl(3,0) and octonions."""

from .quaternion import Quaternion
from .cl3 import Cl3Element
from .products import (
    clifford_product,
    octonion_twin_product,
    associator,
    commutator,
)
from .basis import (
    one,
    qi,
    qj,
    qk,
    ell,
    iell,
    jell,
    kell,
    basis,
)
from .triality import triality_form, inner_product, real_part

__all__ = [
    "Quaternion",
    "Cl3Element",
    "clifford_product",
    "octonion_twin_product",
    "associator",
    "commutator",
    "one",
    "qi",
    "qj",
    "qk",
    "ell",
    "iell",
    "jell",
    "kell",
    "basis",
    "triality_form",
    "inner_product",
    "real_part",
]
