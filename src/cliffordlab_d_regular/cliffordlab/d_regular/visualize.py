import numpy as np
import matplotlib.pyplot as plt

__all__ = ["plot_characteristics", "plot_field_abs", "plot_idempotent_components"]

def plot_characteristics(ax=None, lim=2.0, n_lines=9):
    """
    Plot u=0 and v=0 lines and some level sets u=const, v=const.
    """
    if ax is None:
        fig, ax = plt.subplots()
    xs = np.linspace(-lim, lim, 400)
    ax.plot(xs,  xs, lw=1, label="u=0 (y = x)")
    ax.plot(xs, -xs, lw=1, label="v=0 (y = -x)")
    cs = np.linspace(-lim, lim, n_lines)
    for c in cs:
        ax.plot(xs, c - xs, lw=0.5, alpha=0.5)   # u = c => y = c - x
        ax.plot(xs, xs - c, lw=0.5, alpha=0.5)   # v = c => y = x - c
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_aspect("equal")
    ax.legend(loc="upper right")
    return ax

def plot_field_abs(p, q, x, y, ax=None, title="|f|_M pseudo-norm"):
    """
    Plot Minkowski pseudo-norm |f|_M^2 = p^2 - q^2 as an image.
    """
    if ax is None:
        fig, ax = plt.subplots()
    pm = p**2 - q**2
    im = ax.imshow(pm, origin="lower",
                   extent=[x.min(), x.max(), y.min(), y.max()],
                   aspect="auto")
    ax.set_title(title); ax.set_xlabel("x"); ax.set_ylabel("y")
    plt.colorbar(im, ax=ax)
    return ax

def plot_idempotent_components(Fu, Gv, x, y):
    """
    Given component values Fu (on u-grid) and Gv (on v-grid), show them.
    For simplicity, we just imshow on (x,y) extents.
    """
    fig, axs = plt.subplots(1, 2, figsize=(10,4))
    im1 = axs[0].imshow(Fu, origin="lower",
                        extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto")
    axs[0].set_title("F(u) component")
    plt.colorbar(im1, ax=axs[0])
    im2 = axs[1].imshow(Gv, origin="lower",
                        extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto")
    axs[1].set_title("G(v) component")
    plt.colorbar(im2, ax=axs[1])
    for ax in axs:
        ax.set_xlabel("x"); ax.set_ylabel("y")
    return axs
