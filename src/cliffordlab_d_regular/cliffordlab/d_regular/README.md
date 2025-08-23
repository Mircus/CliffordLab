# cliffordlab/d_regular — Cl(1,0) Regularity Toolkit

Minimal, focused tools to explore D-holomorphic regularity for the split-complex (hyperbolic) numbers D = Cl(1,0), with j^2 = +1.

## Features
- **Builders**: create D-holomorphic fields via the idempotent form f = F(u)e_+ + G(v)e_-.
- **Regularity checks**: numeric CR residuals, Morera-style rectangle integrals, identity tests.
- **Kernels**: principal-value style approximations for (u-α)^{-1}, (v-β)^{-1}.
- **Visualization**: characteristic lines (u=0, v=0), pseudo-norm heatmaps.

## Install (editable)
```bash
pip install -e .
```

*(This is a pure-python subpackage—no compiled deps.)*

## Quickstart
```python
import numpy as np
from cliffordlab.d_regular.builders import from_FG
from cliffordlab.d_regular.regularity import is_d_holomorphic

# grid
x = np.linspace(-2, 2, 201); y = np.linspace(-2, 2, 201)
X, Y = np.meshgrid(x, y)

# F(u), G(v)
F = lambda u: np.exp(0.2*u)
G = lambda v: np.cos(v)

# build and test
p, q = from_FG(F, G, X, Y)
ok, resid = is_d_holomorphic(p, q, X, Y, tol=1e-2)
print("D-holomorphic?", ok, "residual:", resid)
```

## Notes
- We represent D-values by real pairs (p, q) with f = p + j q.
- Idempotent split: F = p + q, G = p - q; combine: p=(F+G)/2, q=(F-G)/2.
- Zero divisors live on u=0 or v=0; avoid them for inverses/logs.

## License
MIT
