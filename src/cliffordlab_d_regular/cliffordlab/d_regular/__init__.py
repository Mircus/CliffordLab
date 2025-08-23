"""
cliffordlab.d_regular
=====================

Tools for exploring regularity and D-holomorphicity in Cl(1,0) (split-complex / hyperbolic numbers).
This package is *numeric-first*: it provides checks, builders, and visualization helpers to explore
the calculus f = F(u) e_+ + G(v) e_- with u = x + y, v = x - y, and j^2 = +1.

Conventions:
- A D-valued function f(x,y) is represented by two real components (p, q) with f = p + j q.
- Idempotents: e_+ = (1 + j)/2, e_- = (1 - j)/2.
- u = x + y, v = x - y.

Modules:
- hyperbolic: basic algebra, utilities, special functions
- builders: construct D-holomorphic functions from one-variable F, G
- regularity: numeric CR checks, Morera rectangles, identity tests
- cauchy_kernel: principal-value regularizations for 1/(u-a) and 1/(v-b)
- visualize: plotting helpers for characteristics and fields
"""
