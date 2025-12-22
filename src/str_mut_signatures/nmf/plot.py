from __future__ import annotations

from typing import Any, Iterable, Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

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

    if n_components < 2:
        raise ValueError("n_components must be at least 2 to plot PC1 vs PC2.")

    X = matrix.to_numpy(dtype=float)

    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X)

    cols = [f"PC{i+1}" for i in range(n_components)]
    coords = pd.DataFrame(X_pca, index=matrix.index, columns=cols)

    return coords, pca.explained_variance_ratio_


def cluster_rows(
    matrix: pd.DataFrame,
    max_clusters: int = 6,
    random_state: int | None = 0,
) -> np.ndarray | None:
    """
    Cluster rows of a matrix using KMeans and select best k via silhouette score.

    Parameters
    ----------
    matrix : pandas.DataFrame
        Rows are samples to cluster (e.g. PCA coordinates).

    max_clusters : int, default 6
        Maximum number of clusters to consider.

    random_state : int or None, default 0
        Random seed for KMeans.

    Returns
    -------
    labels : np.ndarray or None
        Cluster labels (0..k-1) in the same order as matrix.index,
        or None if clustering is not possible.
    """
    n_samples = matrix.shape[0]
    if n_samples < 3:  # not enough for silhouette
        return None

    X = matrix.to_numpy()
    best_k = None
    best_score = -np.inf
    best_labels = None

    for k in range(2, min(max_clusters, n_samples - 1) + 1):
        km = KMeans(n_clusters=k, n_init="auto", random_state=random_state)
        labels = km.fit_predict(X)
        if len(np.unique(labels)) < 2:
            continue
        try:
            score = silhouette_score(X, labels)
        except ValueError:
            continue
        if score > best_score:
            best_score = score
            best_k = k
            best_labels = labels

    return best_labels if best_k is not None else None


def plot_pca_samples(
    result: NMFResult,
    *,
    matrix: pd.DataFrame | None = None,
    n_components: int = 2,
    ax: plt.Axes | None = None,
    title: str | None = None,
    alpha: float = 0.8,
    s: float = 30.0
) -> tuple[pd.DataFrame, np.ndarray, plt.Axes]:
    """
    Run PCA on an NMF result (typically exposures) and plot PC1 vs PC2.

    This is the main entry point:
    - extracts a samples x features matrix from `result`
      (by default `result.exposures`),
    - computes PCA,
    - color samples by NMFResult.groups
    - returns PCA coordinates, variance explained, cluster labels and axes.

    Parameters
    ----------
    result : NMFResult
        Output of `run_nmf`. Must have `.exposures` as a pandas DataFrame,
        unless a `matrix` is explicitly passed.

    matrix : pandas.DataFrame or None, default None
        Optional matrix to use instead of `result.exposures`.
        Must be samples x features. If None, uses `result.exposures`.

    n_components : int, default 2
        Number of principal components to compute. Must be >= 2.

    ax : matplotlib.axes.Axes or None, default None
        Existing axes to plot on. If None, a new figure/axes is created.

    title : str or None, default None
        Plot title. If None, a default title is generated.

    alpha : float, default 0.8
        Point transparency.

    s : float, default 30.0
        Point size.

    Returns
    -------
    coords : pandas.DataFrame
        PCA coordinates (PC1, PC2, ...).

    explained_variance_ratio_ : np.ndarray
        Fraction of variance explained by each component.

    ax : matplotlib.axes.Axes
        Axes with the PCA scatter plot.
    """
    # choose matrix
    if matrix is None:
        if not hasattr(result, "exposures"):
            raise AttributeError(
                "result has no attribute 'exposures'. "
                "Pass a `matrix` argument explicitly if needed."
            )
        matrix = result.exposures

    if not isinstance(matrix, pd.DataFrame):
        raise TypeError("matrix must be a pandas.DataFrame.")

    # run PCA
    coords, var_ratio = compute_pca(matrix, n_components=n_components)

    if "PC1" not in coords.columns or "PC2" not in coords.columns:
        raise ValueError("coords must contain 'PC1' and 'PC2' columns.")

    groups_df = result.groups
    if "group" not in groups_df.columns:
        raise ValueError("result.groups must contain a column named 'group'.")

    # align
    plot_df = coords.join(groups_df[["group"]], how="inner")
    if plot_df.empty:
        raise ValueError("No samples left to plot after aligning PCA coords and groups.")

    # prepare axes
    if ax is None:
        fig, ax = plt.subplots()

    # plot per group to get a correct categorical legend
    for g, sub in plot_df.groupby("group", dropna=False):
        ax.scatter(
            sub["PC1"].to_numpy(),
            sub["PC2"].to_numpy(),
            alpha=alpha,
            s=s,
            label=str(g),
        )

    ax.legend(title="Group", fontsize="small")

    # axis labels with % variance
    pc1_var = var_ratio[0] * 100 if len(var_ratio) > 0 else None
    pc2_var = var_ratio[1] * 100 if len(var_ratio) > 1 else None

    if pc1_var is not None:
        ax.set_xlabel(f"PC1 ({pc1_var:.1f}%)")
    else:
        ax.set_xlabel("PC1")

    if pc2_var is not None:
        ax.set_ylabel(f"PC2 ({pc2_var:.1f}%)")
    else:
        ax.set_ylabel("PC2")

    if title is None:
        title = "PCA of NMF exposures"
    ax.set_title(title)

    return coords, var_ratio, ax


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

def order_by_group_and_total(
    exposures: pd.DataFrame,
    groups: pd.DataFrame
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Align exposures with groups, then reorder rows by group and within each group
    by total exposure (descending).

    Returns
    -------
    ordered_exposures : pd.DataFrame
        Exposures aligned + sorted.

    ordered_groups : np.ndarray
        Group labels in the same order as ordered_exposures.index.
    """
    if not isinstance(exposures, pd.DataFrame):
        raise TypeError("exposures must be a pandas.DataFrame.")
    if not isinstance(groups, pd.DataFrame):
        raise TypeError("groups must be a pandas.DataFrame.")
    if "group" not in groups.columns:
        raise ValueError("groups must contain column 'group'.")

    # Align
    tmp = exposures.join(groups[["group"]], how="inner")
    if tmp.empty:
        raise ValueError("No overlap between exposures.index and groups.index.")

    # Compute totals on exposures only
    totals = tmp.drop(columns=["group"]).sum(axis=1)

    # Stable sort:
    # 1) by group
    # 2) by total exposure (desc)
    tmp["_total"] = totals
    tmp = tmp.sort_values(by=["group", "_total"], ascending=[True, False], kind="mergesort")

    ordered_groups = tmp["group"].to_numpy()
    ordered_exposures = tmp.drop(columns=["group", "_total"])
    return ordered_exposures, ordered_groups


def chunk_indices(n: int, max_per_chunk: int | None) -> list[np.ndarray]:
    if max_per_chunk is None or n <= max_per_chunk:
        return [np.arange(n)]
    return [np.arange(start, min(n, start + max_per_chunk)) for start in range(0, n, max_per_chunk)]


def plot_single_exposure_panel(
    exposures: pd.DataFrame,
    group_labels: np.ndarray | None = None,
    stacked: bool = True,
    figsize: tuple[float, float] | None = None,
    title: str = "Sample exposures",
) -> plt.Figure:
    """
    Plot one panel (subset of samples) of exposures.
    If group_labels is provided, insert gaps between groups.
    """
    n_samples, n_sig = exposures.shape

    if figsize is None:
        figsize = (max(6.0, 0.4 * n_samples), 4.0)

    fig, ax = plt.subplots(figsize=figsize)

    # x positions with optional gaps between groups
    if group_labels is None or n_samples == 0:
        x = np.arange(n_samples, dtype=float)
    else:
        x = np.zeros(n_samples, dtype=float)
        x[0] = 0.0
        gap = 1.0
        for i in range(1, n_samples):
            x[i] = x[i - 1] + 1.0
            if group_labels[i] != group_labels[i - 1]:
                x[i] += gap

    sig_cols = list(exposures.columns)

    if stacked:
        bottom = np.zeros(n_samples)
        for col in sig_cols:
            vals = exposures[col].to_numpy()
            ax.bar(x, vals, bottom=bottom, label=str(col))
            bottom += vals
    else:
        width = 0.8 / max(n_sig, 1)
        for i, col in enumerate(sig_cols):
            vals = exposures[col].to_numpy()
            ax.bar(x + i * width, vals, width=width, label=str(col))
        ax.set_xlim(x.min() - 0.5, x.max() + 0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(exposures.index, rotation=90, fontsize="x-small")

    ax.set_ylabel("Exposure")
    ax.legend(fontsize="small", title="Signatures")
    ax.set_title(title)

    # Add group labels on top (optional, but useful)
    if group_labels is not None and n_samples > 0:
        # annotate group at the start of each group block
        starts = [0] + [i for i in range(1, n_samples) if group_labels[i] != group_labels[i - 1]]
        for i in starts:
            ax.text(x[i], ax.get_ylim()[1], str(group_labels[i]),
                    ha="center", va="bottom", fontsize="x-small", rotation=0)

    fig.tight_layout()
    return fig


def plot_exposures(
    result: NMFResult,
    *,
    stacked: bool = True,
    figsize: tuple[float, float] | None = None,
    max_samples_per_fig: int | None = None,
    plot: Literal["both", "absolute", "proportion"] = "both",
) -> dict[str, list[plt.Figure]]:
    """
    Plot sample exposures to signatures, ordered/grouped by `result.groups["group"]`
    (no new clustering).

    Notes
    -----
    - Always aligns exposures with groups first (inner join).
    - Inserts gaps between groups in the x-axis.
    """
    exp = result.exposures
    if not isinstance(exp, pd.DataFrame):
        raise TypeError("result.exposures must be a pandas.DataFrame.")

    groups = result.groups
    if not isinstance(groups, pd.DataFrame) or groups.empty:
        # fallback: single group for all samples
        groups = pd.DataFrame({"group": "1"}, index=exp.index)
    elif "group" not in groups.columns:
        raise ValueError(f"result.groups must contain column {'group'!r}.")

    # ---- align + order by group and total ----
    exp_ord, group_labels = order_by_group_and_total(exp, groups)

    # proportional exposures (avoid division by zero)
    totals = exp_ord.sum(axis=1)
    prop = exp_ord.div(totals.replace(0, np.nan), axis=0).fillna(0.0)

    figs: dict[str, list[plt.Figure]] = {}
    n_samples = exp_ord.shape[0]
    chunks = chunk_indices(n_samples, max_samples_per_fig)

    if plot in ("both", "absolute"):
        abs_figs: list[plt.Figure] = []
        for chunk_idx in chunks:
            exp_chunk = exp_ord.iloc[chunk_idx]
            gl_chunk = group_labels[chunk_idx] if group_labels is not None else None
            abs_figs.append(
                plot_single_exposure_panel(
                    exp_chunk,
                    group_labels=gl_chunk,
                    stacked=stacked,
                    figsize=figsize,
                    title="Sample exposures (absolute)",
                )
            )
        figs["absolute"] = abs_figs

    if plot in ("both", "proportion"):
        prop_figs: list[plt.Figure] = []
        for chunk_idx in chunks:
            prop_chunk = prop.iloc[chunk_idx]
            gl_chunk = group_labels[chunk_idx] if group_labels is not None else None
            fig = plot_single_exposure_panel(
                prop_chunk,
                group_labels=gl_chunk,
                stacked=True,
                figsize=figsize,
                title="Sample exposures (proportions)",
            )
            ax = fig.axes[0]
            ax.set_ylabel("Proportion")
            ax.set_ylim(0, 1)
            prop_figs.append(fig)
        figs["proportion"] = prop_figs

    return figs
