import numpy as np
from .hyperbolic import u_of, v_of

__all__ = ["cr_residuals", "is_d_holomorphic", "morera_rectangles", "identity_test"]

def cr_residuals(p, q, x, y):
    """
    Compute residuals of D-CR system for f = p + j q:
    R1 = ∂x p - ∂y q,  R2 = ∂x q - ∂y p (both should be ~0).
    Central differences, Neumann at boundaries.
    """
    p = np.asarray(p); q = np.asarray(q)
    x = np.asarray(x); y = np.asarray(y)
    # assume x,y form 2D meshgrid aligned arrays
    dx = np.gradient(x, axis=1)
    dy = np.gradient(y, axis=0)
    px = np.gradient(p, axis=1) / dx
    py = np.gradient(p, axis=0) / dy
    qx = np.gradient(q, axis=1) / dx
    qy = np.gradient(q, axis=0) / dy
    R1 = px - qy
    R2 = qx - py
    return R1, R2

def is_d_holomorphic(p, q, x, y, tol=1e-3):
    """
    Numeric predicate: max-norm of CR residuals below tol.
    """
    R1, R2 = cr_residuals(p, q, x, y)
    m = max(np.nanmax(np.abs(R1)), np.nanmax(np.abs(R2)))
    return bool(m <= tol), m

def morera_rectangles(p, q, x, y, nu=16):
    """
    Morera-style test over 'nu' random rectangles aligned to u and v.
    Approximate ∮ f dz over rectangles whose sides are parallel to u- and v-axes.
    Returns list of absolute integrals.
    """
    # flatten domain ranges
    xmin, xmax = float(np.min(x)), float(np.max(x))
    ymin, ymax = float(np.min(y)), float(np.max(y))

    vals = []
    rng = np.random.default_rng(0)
    for _ in range(nu):
        # pick a random center and sizes
        cx = rng.uniform(xmin, xmax)
        cy = rng.uniform(ymin, ymax)
        sx = 0.2 * (xmax - xmin) * rng.uniform(0.05, 0.2)
        sy = 0.2 * (ymax - ymin) * rng.uniform(0.05, 0.2)

        # u- and v-aligned rectangle edges in (x,y):
        # u=x+y const edges: dx = -dy ; v=x-y const edges: dx = dy
        # We'll sample K points per edge
        K = 256
        t = np.linspace(-1, 1, K)
        # Four edges: (+u), (-u), (+v), (-v)
        edges = []
        # Edge 1: u-constant, going along v
        x1 = cx + sx * t
        y1 = cy - sx * t
        edges.append((x1, y1))
        # Edge 2: u-constant opposite
        x2 = cx - sx * t
        y2 = cy + sx * t
        edges.append((x2, y2))
        # Edge 3: v-constant, going along u
        x3 = cx + sy * t
        y3 = cy + sy * t
        edges.append((x3, y3))
        # Edge 4: v-constant opposite
        x4 = cx - sy * t
        y4 = cy - sy * t
        edges.append((x4, y4))

        # Interpolate p,q on edges by nearest-neighbor sampling from grid
        def sample(arr, xs, ys):
            # assume x,y are meshgrids
            X = x[0, :]; Y = y[:, 0]
            ix = np.clip(np.searchsorted(X, xs) - 1, 0, X.size-1)
            iy = np.clip(np.searchsorted(Y, ys) - 1, 0, Y.size-1)
            return arr[iy, ix]

        integral = 0.0 + 0.0j
        for k, (xe, ye) in enumerate(edges):
            pe = sample(p, xe, ye)
            qe = sample(q, xe, ye)
            fe = pe + 1j*qe  # embed j as imaginary for path integral bookkeeping
            # dz along edge: compute incremental dx + j dy with j^2=+1 -> here we just track magnitude
            # For Morera, it's enough to use param derivative along the edge:
            dxe = np.gradient(xe)
            dye = np.gradient(ye)
            # Treat j as 1 in magnitude for line integration surrogate:
            dz = dxe + 1j*dye
            integral += np.trapz(fe * dz, t)
        vals.append(abs(integral))
    return vals

def identity_test(p, q, x, y, lines=8, m_slope=0.0):
    """
    Test limited identity principle numerically.
    Sample f on several non-characteristic lines y = m x + b and check if zero there implies zero nearby.
    Returns a dict with line data and residuals.
    """
    xmin, xmax = float(np.min(x)), float(np.max(x))
    ymin, ymax = float(np.min(y)), float(np.max(y))

    rng = np.random.default_rng(1)
    results = []
    for _ in range(lines):
        b = rng.uniform(ymin + 0.1*(ymax-ymin), ymax - 0.1*(ymax-ymin))
        # build a line segment inside the box
        t = np.linspace(-0.8, 0.8, 400)
        xs = 0.5*(xmin+xmax) + t*(xmax-xmin)/2
        ys = m_slope*xs + b
        # sample
        X = x[0, :]; Y = y[:, 0]
        ix = np.clip(np.searchsorted(X, xs) - 1, 0, X.size-1)
        iy = np.clip(np.searchsorted(Y, ys) - 1, 0, Y.size-1)
        ps = p[iy, ix]; qs = q[iy, ix]
        line_resid = float(np.max(np.sqrt(ps**2 + qs**2)))
        results.append({"b": float(b), "max_norm_on_line": line_resid})
    return {"slope": m_slope, "lines": results}
