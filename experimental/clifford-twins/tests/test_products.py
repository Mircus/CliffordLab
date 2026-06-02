from clifford_twins import (
    clifford_product,
    octonion_twin_product,
    associator,
)
from clifford_twins.basis import one, qi, qj, qk, ell
from clifford_twins.identities import (
    is_associative_on_basis,
    is_alternative_on_basis,
)


def test_clifford_product_associative_on_basis():
    assert is_associative_on_basis(clifford_product)


def test_octonion_twin_not_associative_on_basis():
    assert not is_associative_on_basis(octonion_twin_product)


def test_octonion_twin_alternative_on_basis():
    assert is_alternative_on_basis(octonion_twin_product)


def test_active_volume_crossing_rule():
    # ell * i = conjugate(i) * ell = -i * ell
    lhs = octonion_twin_product(ell, qi)
    rhs = -octonion_twin_product(qi, ell)
    assert lhs.approx_eq(rhs)


def test_nonassociativity_witness():
    lhs = octonion_twin_product(octonion_twin_product(ell, qi), qj)
    rhs = octonion_twin_product(ell, octonion_twin_product(qi, qj))

    assert not lhs.approx_eq(rhs)
    assert not associator(octonion_twin_product, ell, qi, qj).is_zero()


def test_norm_multiplicativity_for_basis_sample():
    x = qi + ell
    y = qj + qk
    xy = octonion_twin_product(x, y)
    assert abs(xy.norm2() - x.norm2() * y.norm2()) < 1e-9
