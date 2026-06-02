import random

from clifford_twins import (
    octonion_twin_product,
    clifford_product,
)
from clifford_twins.basis import basis, one
from clifford_twins.cl3 import Cl3Element


def rand_quat(scale=3.0):
    import random
    return (random.uniform(-scale, scale),
            random.uniform(-scale, scale),
            random.uniform(-scale, scale),
            random.uniform(-scale, scale))


def rand_elem(scale=3.0):
    ar, ax, ay, az = rand_quat(scale)
    br, bx, by, bz = rand_quat(scale)
    from clifford_twins.quaternion import Quaternion
    return Cl3Element(Quaternion(ar, ax, ay, az), Quaternion(br, bx, by, bz))


def test_octonion_conjugation_involution_and_anti_automorphism():
    for _ in range(50):
        x = rand_elem()
        y = rand_elem()
        xx = x.conj_octonion().conj_octonion()
        assert x.approx_eq(xx)
        # anti-automorphism: \overline{xy} = \bar y \bar x
        xy = octonion_twin_product(x, y)
        lhs = xy.conj_octonion()
        rhs = octonion_twin_product(y.conj_octonion(), x.conj_octonion())
        assert lhs.approx_eq(rhs, tol=1e-7)


def test_octonion_quadratic_form_and_norm_multiplicativity_random():
    for _ in range(50):
        x = rand_elem()
        y = rand_elem()
        # x * \bar x = \bar x * x = ||x||^2
        xxb = octonion_twin_product(x, x.conj_octonion())
        bbx = octonion_twin_product(x.conj_octonion(), x)
        n2 = x.norm2()
        assert abs(xxb.a.r - n2) < 1e-7 and xxb.b.approx_eq(xxb.b.__class__(), tol=1e-9)
        assert abs(bbx.a.r - n2) < 1e-7 and bbx.b.approx_eq(bbx.b.__class__(), tol=1e-9)
        # multiplicativity
        xy = octonion_twin_product(x, y)
        assert abs(xy.norm2() - x.norm2() * y.norm2()) < 1e-6


def test_associator_real_part_is_zero():
    # Re((xy)z) == Re(x(yz)) should hold in octonions
    for _ in range(50):
        x = rand_elem()
        y = rand_elem()
        z = rand_elem()
        left = octonion_twin_product(octonion_twin_product(x, y), z).real()
        right = octonion_twin_product(x, octonion_twin_product(y, z)).real()
        assert abs(left - right) < 1e-7


def test_alternativity_random():
    for _ in range(50):
        x = rand_elem()
        y = rand_elem()
        a1 = octonion_twin_product(octonion_twin_product(x, x), y)
        a2 = octonion_twin_product(x, octonion_twin_product(x, y))
        b1 = octonion_twin_product(y, octonion_twin_product(x, x))
        b2 = octonion_twin_product(octonion_twin_product(y, x), x)
        assert a1.approx_eq(a2, tol=1e-7)
        assert b1.approx_eq(b2, tol=1e-7)


def test_flexibility_random():
    for _ in range(50):
        x = rand_elem()
        y = rand_elem()
        lhs = octonion_twin_product(octonion_twin_product(x, y), x)
        rhs = octonion_twin_product(x, octonion_twin_product(y, x))
        assert lhs.approx_eq(rhs, tol=1e-7)


def test_clifford_vs_octonion_on_even_core():
    # When b=0, star reduces to H product and matches Clifford product too
    from clifford_twins.quaternion import Quaternion
    import random
    for _ in range(20):
        a = Cl3Element(Quaternion(*rand_quat()), Quaternion())
        c = Cl3Element(Quaternion(*rand_quat()), Quaternion())
        star = octonion_twin_product(a, c)
        cliff = clifford_product(a, c)
        assert star.approx_eq(cliff, tol=1e-9)


def test_square_of_pure_ell_component_is_negative_real():
    from clifford_twins.quaternion import Quaternion
    for q in [
        Quaternion(0, 0, 0, 0),
        Quaternion(1, 0, 0, 0),
        Quaternion(0, 2, -3, 1),
    ]:
        x = Cl3Element(Quaternion(), q)
        xx = octonion_twin_product(x, x)
        # equals -(\bar q q) which is -||q||^2
        expected = - (q.conj() * q).r
        assert abs(xx.a.r - expected) < 1e-9 and xx.b.approx_eq(Quaternion(), tol=1e-9)
