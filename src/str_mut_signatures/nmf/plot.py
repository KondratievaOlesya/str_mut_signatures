from __future__ import annotations

from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from .nmf import NMFResult

# ---------------------------------------------------------------------
# PCA helpers
# ---------------------------------------------------------------------


def compute_pca(
    matrix: pd.DataFrame,
    n_components: int = 2,
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Compute PCA on a samples x features matrix (e.g. exposures).

    Parameters
    ----------
    matrix : pandas.DataFrame
        Numeric matrix:
        - rows   : samples
        - columns: features or signatures.

    n_components : int, default 2
        Number of principal components to compute.

    Returns
    -------
    coords : pandas.DataFrame
        PCA coordinates with:
        - index  : same as matrix.index
        - columns: PC1, PC2, ..., PC{n_components}

    explained_variance_ratio_ : np.ndarray
        1D array of length n_components with the fraction of variance
        explained by each component.
    """
    if not isinstance(matrix, pd.DataFrame):
        raise TypeError("matrix must be a pandas.DataFrame.")

    if matrix.empty:
        raise ValueError("matrix is empty; cannot compute PCA.")

    if not all(np.issubdtype(dtype, np.number) for dtype in matrix.dtypes.values):
        raise TypeError("matrix must contain only numeric values for PCA.")

    X = matrix.to_numpy(dtype=float)

    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X)

    cols = [f"PC{i+1}" for i in range(n_components)]
    coords = pd.DataFrame(X_pca, index=matrix.index, columns=cols)

    return coords, pca.explained_variance_ratio_


def plot_pca_samples(
    coords: pd.DataFrame,
    labels: Iterable[Any] | None = None,
    ax: plt.Axes | None = None,
    title: str | None = None,
    alpha: float = 0.8,
    s: float = 30.0,
) -> plt.Axes:
    """
    Scatter plot of samples in PCA space (PC1 vs PC2).

    Parameters
    ----------
    coords : pandas.DataFrame
        Output of `compute_pca` (at least PC1, PC2 columns).

    labels : iterable or None, default None
        Optional labels for coloring points. If provided, must have
        the same length as coords.index. Values are used as categories.

    ax : matplotlib.axes.Axes or None, default None
        Existing axes to plot on. If None, a new figure/axes is created.

    title : str or None, default None
        Plot title.

    alpha : float, default 0.8
        Point transparency.

    s : float, default 30.0
        Point size.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes with the plot.
    """
    if ax is None:
        fig, ax = plt.subplots()

    if "PC1" not in coords.columns or "PC2" not in coords.columns:
        raise ValueError("coords must contain 'PC1' and 'PC2' columns.")

    x = coords["PC1"].to_numpy()
    y = coords["PC2"].to_numpy()

    if labels is None:
        ax.scatter(x, y, alpha=alpha, s=s)
    else:
        labels = list(labels)
        if len(labels) != len(coords):
            raise ValueError("labels must have the same length as coords.index.")

        labels_arr = np.array(labels)
        unique = np.unique(labels_arr)

        for lab in unique:
            mask = labels_arr == lab
            ax.scatter(
                x[mask],
                y[mask],
                alpha=alpha,
                s=s,
                label=str(lab),
            )
        ax.legend(title="Group", fontsize="small")

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    if title is not None:
        ax.set_title(title)

    return ax


# ---------------------------------------------------------------------
# Signature plots
# ---------------------------------------------------------------------


def plot_signatures(
    result: NMFResult,
    top_n: int = 20,
    signatures: list[int] | list[str] | None = None,
    figsize: tuple[float, float] | None = None,
    sharey: bool = False,
) -> plt.Figure:
    """
    Plot per-signature barplots of feature loadings.

    Parameters
    ----------
    result : NMFResult
        Output of `run_nmf`.

    top_n : int, default 20
        Number of top features (by loading) to display per signature.

    signatures : list[int] or list[str] or None, default None
        Which signatures to plot.
        - If None, plot all.
        - If list[int], interpreted as 1-based indices (1..K).
        - If list[str], must match column names in result.signatures.

    figsize : (float, float) or None, default None
        Figure size passed to matplotlib. If None, a default is chosen
        based on the number of signatures.

    sharey : bool, default False
        If True, all subplots share the same y-axis.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure.
    """
    sig_df = result.signatures

    # Determine which signatures to plot
    if signatures is None:
        sig_cols = list(sig_df.columns)
    elif isinstance(signatures[0], int):
        # 1-based indices
        sig_cols = [sig_df.columns[i - 1] for i in signatures]
    else:
        sig_cols = list(signatures)

    n_sig = len(sig_cols)
    if n_sig == 0:
        raise ValueError("No signatures selected to plot.")

    # Figure layout
    if figsize is None:
        figsize = (4 * n_sig, 4)

    fig, axes = plt.subplots(1, n_sig, figsize=figsize, sharey=sharey)
    if n_sig == 1:
        axes = [axes]

    for ax, col in zip(axes, sig_cols):
        s = sig_df[col]

        # pick top_n features by loading
        s_sorted = s.sort_values(ascending=False).head(top_n)

        ax.bar(range(len(s_sorted)), s_sorted.values)
        ax.set_xticks(range(len(s_sorted)))
        ax.set_xticklabels(s_sorted.index, rotation=90, fontsize="x-small")
        ax.set_title(str(col))
        ax.set_ylabel("Loading")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------
# Exposure / sample plots
# ---------------------------------------------------------------------


def plot_exposures(
    result: NMFResult,
    sort_by: int | None = None,
    stacked: bool = True,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """
    Plot sample exposures to signatures.

    Parameters
    ----------
    result : NMFResult
        Output of `run_nmf`.

    sort_by : int or None, default None
        If not None, sort samples by exposure to this signature.
        - Interpreted as 1-based index (1..K).

    stacked : bool, default True
        If True, draw stacked barplot of exposures per sample.
        If False, draw grouped bars (one bar per signature side-by-side).

    figsize : (float, float) or None, default None
        Figure size. If None, choose a basic default depending on
        number of samples.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure.
    """
    exp = result.exposures.copy()

    if sort_by is not None:
        n_sig = exp.shape[1]
        if not (1 <= sort_by <= n_sig):
            raise ValueError(f"sort_by must be between 1 and {n_sig} (1-based).")
        col = exp.columns[sort_by - 1]
        exp = exp.sort_values(by=col, ascending=False)

    n_samples, n_sig = exp.shape
    if figsize is None:
        figsize = (max(6.0, 0.4 * n_samples), 4.0)

    fig, ax = plt.subplots(figsize=figsize)

    x = np.arange(n_samples)
    sig_cols = list(exp.columns)

    if stacked:
        bottom = np.zeros(n_samples)
        for col in sig_cols:
            vals = exp[col].to_numpy()
            ax.bar(x, vals, bottom=bottom, label=str(col))
            bottom += vals
    else:
        # grouped bars
        width = 0.8 / n_sig
        for i, col in enumerate(sig_cols):
            vals = exp[col].to_numpy()
            ax.bar(x + i * width, vals, width=width, label=str(col))
        ax.set_xlim(-0.2, n_samples - 0.2 + 0.8)

    ax.set_xticks(x + (0 if stacked else 0.4))
    ax.set_xticklabels(exp.index, rotation=90, fontsize="x-small")

    ax.set_ylabel("Exposure")
    ax.legend(fontsize="small", title="Signatures")
    ax.set_title("Sample exposures")

    fig.tight_layout()
    return fig
