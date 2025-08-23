"""
Demo: build a D-holomorphic function, check CR, run Morera rectangles,
and visualize characteristics.
"""
import numpy as np
import matplotlib.pyplot as plt

from cliffordlab.d_regular.builders import from_FG, sqrt_abs_regularized
from cliffordlab.d_regular.regularity import cr_residuals, is_d_holomorphic, morera_rectangles
from cliffordlab.d_regular.visualize import plot_characteristics, plot_field_abs

# domain grid
nx = ny = 201
x = np.linspace(-2, 2, nx)
y = np.linspace(-2, 2, ny)
X, Y = np.meshgrid(x, y)

# Build F,G
F = sqrt_abs_regularized(eps=1e-4)
G = lambda t: np.exp(0.5 * t)

# Construct f = F(u)e_+ + G(v)e_-
p, q = from_FG(F, G, X, Y)

# Check CR residuals
R1, R2 = cr_residuals(p, q, X, Y)
flag, maxres = is_d_holomorphic(p, q, X, Y, tol=1e-2)
print("D-holomorphic (numeric)?", flag, "max residual:", maxres)

# Morera rectangles
vals = morera_rectangles(p, q, X, Y, nu=12)
print("Morera rectangle integrals (abs):", np.round(vals, 4))
print("Median:", float(np.median(vals)))

# Visualize
ax = plot_characteristics(lim=2.0)
plot_field_abs(p, q, X, Y, ax=ax, title="|f|_M^2 for demo f")
plt.tight_layout()
plt.show()
