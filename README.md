
 <img src="./cliffordlab.svg" width="650" alt="CliffordLab logo" />


# CliffordLab

Minimal, practical tools for experimenting with **unified slice regularity** on real Clifford algebras \(Cl(p,q)\).

- `cliffordlab.py` – tiny, self-contained GA core + unified-slice helpers (no deps)
- `cliffordlab_sym.py` – **SymPy** exact CR checks and slice exponentials
- `cliffordlab_pygae.py` – helpers built on **pygae/clifford** (numeric GA)


## Install

```bash
# core lab has no deps
python -m pip install sympy clifford  # optional backends
# CliffordLab
