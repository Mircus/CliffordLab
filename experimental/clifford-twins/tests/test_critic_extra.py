import itertools as it
import math

from clifford_twins import octonion_twin_product
from clifford_twins.basis import basis
from clifford_twins.triality import real_part
from clifford_twins.identities import associator


B = basis()


def mul(x, y):
    return octonion_twin_product(x, y)


def test_moufang_identity_basis_triples():
    # (xy)(zx) = x(yz)x for octonions (one of Moufang identities)
    for x, y, z in it.product(B, repeat=3):
        lhs = mul(mul(x, y), mul(z, x))
        rhs = mul(x, mul(mul(y, z), x))
        assert lhs.approx_eq(rhs), (x, y, z, lhs, rhs)


def test_real_part_of_associator_is_zero_on_random_samples():
    # For octonions, Re((xy)z) = Re(x(yz))
    # Test over all basis triples suffices by trilinearity of Re and bilinearity of product
    for x, y, z in it.product(B, repeat=3):
        lhs = real_part(mul(mul(x, y), z))
        rhs = real_part(mul(x, mul(y, z)))
        assert math.isclose(lhs, rhs, abs_tol=1e-9)


def test_triality_form_cyclic_symmetry_sample():
    # t(x,y,z) = Re(x(yz)) is cyclically symmetric in octonions
    for x, y, z in it.product(B, repeat=3):
        t1 = real_part(mul(x, mul(y, z)))
        t2 = real_part(mul(y, mul(z, x)))
        t3 = real_part(mul(z, mul(x, y)))
        assert math.isclose(t1, t2, abs_tol=1e-9)
        assert math.isclose(t2, t3, abs_tol=1e-9)


def test_alternator_is_totally_alternating_on_basis():
    # In alternative algebras, [x,y,z] is alternating: vanishes when two args equal
    for x, y in it.product(B, repeat=2):
        assert associator(octonion_twin_product, x, x, y).is_zero()
        assert associator(octonion_twin_product, x, y, x).is_zero()
        assert associator(octonion_twin_product, y, x, x).is_zero()
