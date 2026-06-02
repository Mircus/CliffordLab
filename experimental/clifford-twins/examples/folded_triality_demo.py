from clifford_twins import octonion_twin_product, triality_form
from clifford_twins.basis import one, qi, qj, qk, ell, iell, jell, kell


def main() -> None:
    print("Folded triality demo")
    print("====================")

    x = qi
    psi = qj
    phi = kell

    print("Triality form t(x,psi,phi)=Re(x*(psi*phi))")
    print("x   = i")
    print("psi = j")
    print("phi = k ell")
    print("t   =", triality_form(x, psi, phi))

    print()
    print("The point of the demo is not to implement Spin(8).")
    print("It only shows the compact octonionic slot:")
    print("  vector face, left-spinor face, right-spinor face")
    print("can all be represented by elements of the same 8D body.")


if __name__ == "__main__":
    main()
