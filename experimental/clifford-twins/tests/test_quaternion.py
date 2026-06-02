from clifford_twins.quaternion import Quaternion


def test_quaternion_units():
    one = Quaternion(1, 0, 0, 0)
    i = Quaternion(0, 1, 0, 0)
    j = Quaternion(0, 0, 1, 0)
    k = Quaternion(0, 0, 0, 1)

    assert (i * i).approx_eq(-one)
    assert (j * j).approx_eq(-one)
    assert (k * k).approx_eq(-one)

    assert (i * j).approx_eq(k)
    assert (j * k).approx_eq(i)
    assert (k * i).approx_eq(j)

    assert (j * i).approx_eq(-k)
