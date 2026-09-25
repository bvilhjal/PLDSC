"""Figures for the multi-ancestry measurements (optional, needs matplotlib).

Three plots carry the live project (STRATEGY §2/§4): the κ heatmap of the
headline table (MULTI_ANCESTRY Table 3), the per-chromosome κ profile the
genome-wide step 1 produces, and the before/after forest plot of published
h² the casualty step 3 produces. They are deliberately thin wrappers: all
numbers are computed with the tested functions in ``pldsc.multiancestry``;
this module only draws arrays.

matplotlib is imported inside the functions, so the core package keeps its
numpy-only dependency; install the extra with ``pip install -e ".[plot]"``.
Every function returns the ``Axes`` it drew on and writes nothing unless
the caller passes ``path``.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike

__all__ = ["kappa_heatmap", "kappa_profile", "forest_plot"]


def _ax(ax):
    import matplotlib.pyplot as plt
    return plt.gca() if ax is None else ax


def _save(ax, path):
    if path is not None:
        ax.figure.savefig(path, dpi=200, bbox_inches="tight")


def kappa_heatmap(
    kappa: ArrayLike,
    row_labels: Sequence[str],
    col_labels: Sequence[str],
    *,
    ax=None,
    title: str | None = None,
    cmap: str = "YlOrRd",
    path: str | None = None,
):
    """Heatmap of κ = ĥ²/h² values (the headline table as a picture).

    ``kappa`` is ``(R, C)``; rows are drawn bottom-up so the first row
    label appears at the bottom. Cells are annotated at three decimals;
    1.0 (unbiased) is marked on the colour bar, since the whole point of
    the figure is how far below it the panel sits.
    """
    K = np.asarray(kappa, dtype=np.float64)
    if K.ndim != 2 or K.shape != (len(row_labels), len(col_labels)):
        raise ValueError("kappa must have shape (len(row_labels), len(col_labels))")
    ax = _ax(ax)
    im = ax.imshow(K, cmap=cmap, vmin=min(K.min(), 0.99), vmax=max(K.max(), 1.0))
    ax.set_xticks(range(K.shape[1]), [str(c) for c in col_labels])
    ax.set_yticks(range(K.shape[0]), [str(r) for r in row_labels])
    for i in range(K.shape[0]):
        for j in range(K.shape[1]):
            ax.text(j, i, f"{K[i, j]:.3f}", ha="center", va="center",
                    color="black" if K[i, j] > (K.min() + 0.5 * (K.max() - K.min())) else "white")
    cb = ax.figure.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label(r"$\kappa = \hat h^2 / h^2$")
    if title:
        ax.set_title(title)
    _save(ax, path)
    return ax


def kappa_profile(
    per_chromosome: ArrayLike,
    chromosomes: Sequence[int],
    *,
    ancestries: Sequence[str] | None = None,
    ax=None,
    title: str | None = None,
    path: str | None = None,
):
    """κ per chromosome, one line per non-reference ancestry (step 1).

    ``per_chromosome`` is ``(n_chromosomes, A)`` for the A non-reference
    ancestries at the rule cell (e.g. EUR 0.75 + X, ρ = 0.8). Draws the
    1.0 reference line and the pre-stated gate (κ ≤ 0.90, STRATEGY §4
    step 1) so the reader can see the gate, not just the lines.
    """
    K = np.asarray(per_chromosome, dtype=np.float64)
    if K.ndim != 2 or K.shape[0] != len(chromosomes):
        raise ValueError("per_chromosome must have one row per chromosome")
    if ancestries is None:
        ancestries = [f"X{i + 1}" for i in range(K.shape[1])]
    if K.shape[1] != len(ancestries):
        raise ValueError("one column per ancestry (or pass ancestries)")
    ax = _ax(ax)
    x = np.asarray(chromosomes)
    for a in range(K.shape[1]):
        ax.plot(x, K[:, a], marker="o", ms=3, lw=1, label=str(ancestries[a]))
    ax.axhline(1.0, color="grey", lw=0.8, ls=":")
    ax.axhline(0.90, color="black", lw=0.8, ls="--")
    ax.text(x[-1], 0.90, " gate κ ≤ 0.90 ", va="center", ha="left", fontsize=8)
    ax.set_xlabel("chromosome")
    ax.set_ylabel(r"$\kappa$ (rule cell: EUR 0.75 + X, $\rho$ = 0.8)")
    ax.legend(frameon=False, ncol=min(4, K.shape[1]))
    if title:
        ax.set_title(title)
    _save(ax, path)
    return ax


def forest_plot(
    labels: Sequence[str],
    published: ArrayLike,
    corrected: ArrayLike,
    se: ArrayLike,
    *,
    ax=None,
    title: str | None = None,
    path: str | None = None,
):
    """Published h² vs the recomputed (meta-score) h², per trait (step 3).

    Each trait gets a row: the published estimate as a hollow marker at
    x=1 and the recomputed estimate with its ±1.96 SE whisker; the line
    between them shows the correction's size and sign. Traits are ordered
    by relative change, largest first, so the figure leads with its
    headline. ``se`` is the SE of the *difference* if you want the
    whiskers to show the move's significance rather than the estimate's
    precision — pass ``se=None`` to omit whiskers.
    """
    pub = np.asarray(published, dtype=np.float64)
    cor = np.asarray(corrected, dtype=np.float64)
    if pub.shape != cor.shape or pub.shape != (len(labels),):
        raise ValueError("published, corrected and labels must match in length")
    order = np.argsort(pub / cor)              # biggest relative loss first
    labels = [labels[i] for i in order]
    pub, cor = pub[order], cor[order]
    ax = _ax(ax)
    y = np.arange(len(labels))[::-1]
    for yi, p, c in zip(y, pub, cor):
        ax.plot([c, p], [yi, yi], color="lightgrey", lw=1, zorder=1)
    ax.scatter(pub, y, marker="s", facecolors="none", edgecolors="black",
               zorder=2, label="published (EUR score)")
    ax.scatter(cor, y, marker="o", color="crimson", zorder=3,
               label="recomputed (meta score)")
    if se is not None:
        s = np.asarray(se, dtype=np.float64)[order]
        ax.hlines(y, cor - 1.96 * s, cor + 1.96 * s, color="crimson", lw=1)
    ax.set_yticks(y, labels)
    ax.set_xlabel(r"$h^2$")
    ax.legend(frameon=False, loc="lower right")
    if title:
        ax.set_title(title)
    _save(ax, path)
    return ax
