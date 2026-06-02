from clifford_twins import (
    clifford_product,
    octonion_twin_product,
    associator,
)
from clifford_twins.basis import qi, qj, ell
from clifford_twins.identities import (
    is_associative_on_basis,
    is_alternative_on_basis,
    associator_table,
)


def main() -> None:
    print("Clifford twins demo")
    print("===================")

    print("Passive Clifford product associative on basis:")
    print(" ", is_associative_on_basis(clifford_product))

    print("Active octonion twin product associative on basis:")
    print(" ", is_associative_on_basis(octonion_twin_product))

    print("Active octonion twin product alternative on basis:")
    print(" ", is_alternative_on_basis(octonion_twin_product))

    lhs = octonion_twin_product(octonion_twin_product(ell, qi), qj)
    rhs = octonion_twin_product(ell, octonion_twin_product(qi, qj))

    print()
    print("Nonassociativity witness:")
    print("  (ell i) j =", lhs)
    print("  ell (i j) =", rhs)
    print("  [ell,i,j] =", associator(octonion_twin_product, ell, qi, qj))

    table = associator_table(octonion_twin_product)
    print()
    print("Number of nonzero basis associators for active product:")
    print(" ", len(table))


if __name__ == "__main__":
    main()
