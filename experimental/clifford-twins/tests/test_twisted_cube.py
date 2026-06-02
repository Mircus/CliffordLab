from clifford_twins.twisted_cube import Twist, dot_mod2, vectors


def test_vectors():
    assert vectors(3) == list(range(8))


def test_dot_mod2():
    assert dot_mod2(0b101, 0b001) == 1
    assert dot_mod2(0b101, 0b010) == 0


def test_trivial_twist_is_associative():
    twist = Twist(3, lambda a, b: 1)
    assert twist.is_associative()
    assert twist.multiply_basis(3, 5) == (1, 6)
