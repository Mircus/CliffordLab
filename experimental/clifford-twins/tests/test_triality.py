from clifford_twins import triality_form
from clifford_twins.basis import one, qi, qj, qk, ell


def test_triality_form_is_trilinear_in_first_slot_sample():
    lhs = triality_form(qi + qj, qk, ell)
    rhs = triality_form(qi, qk, ell) + triality_form(qj, qk, ell)
    assert abs(lhs - rhs) < 1e-9


def test_triality_form_returns_float():
    value = triality_form(qi, qj, qk)
    assert isinstance(value, float)
