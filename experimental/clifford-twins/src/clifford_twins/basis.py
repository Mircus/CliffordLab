from .cl3 import Cl3Element
from .quaternion import Quaternion

# Quaternion basis inside H = Cl^+(3,0)
one = Cl3Element(Quaternion(1, 0, 0, 0))
qi = Cl3Element(Quaternion(0, 1, 0, 0))
qj = Cl3Element(Quaternion(0, 0, 1, 0))
qk = Cl3Element(Quaternion(0, 0, 0, 1))

# Volume slot. Under the passive product this is I.
# Under the active product this is ell.
ell = Cl3Element(Quaternion(), Quaternion(1, 0, 0, 0))

# Remaining basis elements in H ell.
iell = Cl3Element(Quaternion(), Quaternion(0, 1, 0, 0))
jell = Cl3Element(Quaternion(), Quaternion(0, 0, 1, 0))
kell = Cl3Element(Quaternion(), Quaternion(0, 0, 0, 1))


def basis() -> list[Cl3Element]:
    """Standard 8D basis [1,i,j,k,ell,iell,jell,kell]."""
    return [one, qi, qj, qk, ell, iell, jell, kell]
