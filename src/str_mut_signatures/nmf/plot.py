from __future__ import annotations

import re
from typing import Literal

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

    cols = [f"PC{i + 1}" for i in range(n_components)]
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
    cmap: str | None = None,
    s: float = 30.0,
) -> tuple[pd.DataFrame, np.ndarray, plt.Axes]:
    """
    Run PCA on an NMF result (typically exposures) and plot PC1 vs PC2.

    This is the main entry point for PCA visualization:

    - Extract a samples × features matrix from ``result`` (by default
      ``result.exposures``).
    - Compute PCA coordinates.
    - Color samples by ``result.groups``.
    - Plot the first two principal components.

    Parameters
    ----------
    result : NMFResult
        Output of :func:`run_nmf`. Must provide ``.exposures`` as a
        :class:`pandas.DataFrame` unless ``matrix`` is explicitly provided.
    matrix : pandas.DataFrame or None, optional
        Optional matrix to use instead of ``result.exposures``.
        Must be samples × features. If ``None``, uses ``result.exposures``.
    n_components : int, optional
        Number of principal components to compute. Must be >= 2.
        Default is 2.
    ax : matplotlib.axes.Axes or None, optional
        Existing axes to plot on. If ``None``, a new figure and axes are created.
    title : str or None, optional
        Plot title. If ``None``, a default title is generated.
    alpha : float, optional
        Point transparency. Default is 0.8.
    cmap : str or None, optional
        Matplotlib colormap name used for coloring samples.

        - If ``group`` is detected as categorical, the default is ``"tab20"``.
        - If ``group`` is detected as continuous, the default is ``"viridis"``.

        If provided explicitly, this colormap is used in both cases.
    s : float, optional
        Point size. Default is 30.0.

    Returns
    -------
    coords : pandas.DataFrame
        PCA coordinates for each sample (``PC1``, ``PC2``, ...), indexed by sample.
    explained_variance_ratio_ : numpy.ndarray
        Fraction of variance explained by each principal component.
    ax : matplotlib.axes.Axes
        Axes containing the PCA scatter plot.

    Raises
    ------
    ValueError
        If ``n_components < 2`` or if the chosen matrix is empty or non-numeric.
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
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots()
        created_fig = True
    else:
        fig = ax.figure

    values = plot_df["group"]

    # autodetect continuous (numeric) vs categorical
    is_numeric = pd.api.types.is_numeric_dtype(plot_df["group"])
    n_non_null = int(plot_df["group"].notna().sum())

    # use nunique(dropna=True) so NaNs don't inflate uniqueness
    n_unique = int(plot_df["group"].nunique(dropna=True)) if n_non_null > 0 else 0
    unique_frac = (n_unique / n_non_null) if n_non_null > 0 else 0.0

    is_continuous = bool(is_numeric and (n_unique > 10 or unique_frac > 0.30))
    if is_continuous:
        # ---- continuous coloring ----
        cmap_name = cmap if cmap is not None else "viridis"
        cmap_obj = plt.get_cmap(cmap_name)
        vmin = np.nanmin(values.to_numpy(dtype=float))
        vmax = np.nanmax(values.to_numpy(dtype=float))
        norm = plt.Normalize(vmin=vmin, vmax=vmax)

        ax.scatter(
            plot_df["PC1"].to_numpy(),
            plot_df["PC2"].to_numpy(),
            c=cmap_obj(norm(values.to_numpy(dtype=float))),
            alpha=alpha,
            s=s,
        )

        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
        mappable.set_array([])
        fig.colorbar(mappable, ax=ax, label="group")

    else:
        cmap_name = cmap if cmap is not None else "tab20"
        # ---- categorical coloring ----
        uniq = plot_df["group"].dropna().unique()
        n_groups = len(uniq)

        cmap_obj = plt.get_cmap(cmap_name, max(n_groups, 1))
        color_map = {g: cmap_obj(i) for i, g in enumerate(uniq)}

        for g, sub in plot_df.groupby("group", dropna=False):
            color = color_map.get(g, "grey")
            ax.scatter(
                sub["PC1"].to_numpy(),
                sub["PC2"].to_numpy(),
                alpha=alpha,
                s=s,
                color=color,
                label=str(g),
            )

        # legend outside
        ax.legend(
            title="Group",
            fontsize="small",
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
        )
        if created_fig:
            fig.subplots_adjust(right=0.78)

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
    Plot per-signature bar plots of feature loadings.

    This function visualizes the strongest features for each selected
    signature as bar plots, using the signature profiles stored in
    ``result.signatures``.

    Parameters
    ----------
    result : NMFResult
        Output of :func:`run_nmf` containing the signature matrix.
    top_n : int, optional
        Number of top features (by absolute loading) to display per signature.
        Default is 20.
    signatures : list[int] or list[str] or None, optional
        Which signatures to plot.

        - If ``None``, all signatures are plotted.
        - If ``list[int]``, values are interpreted as 1-based indices
          (``1 .. K``).
        - If ``list[str]``, values must match column names in
          ``result.signatures``.
    figsize : tuple[float, float] or None, optional
        Figure size passed to matplotlib. If ``None``, a default size is chosen
        based on the number of signatures.
    sharey : bool, optional
        If ``True``, all subplots share the same y-axis. Default is ``False``.

    Returns
    -------
    matplotlib.figure.Figure
        The created matplotlib figure containing the signature bar plots.

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


def build_exposure_table(
    exposures: pd.DataFrame, groups: pd.DataFrame | None = None
) -> pd.DataFrame:
    """
    Return a single dataframe with exposures + a 'group' column (aligned by index).

    - Inner-joins on sample index to guarantee exposure rows match group rows.
    - If groups is missing/empty, assigns all samples to default_group.
    """
    if not isinstance(exposures, pd.DataFrame):
        raise TypeError("exposures must be a pandas.DataFrame.")

    exp = exposures.copy()

    if groups is None or not isinstance(groups, pd.DataFrame) or groups.empty:
        grp = pd.DataFrame({"group": "1"}, index=exp.index)
    else:
        if "group" not in groups.columns:
            raise ValueError('groups must contain column "group".')
        grp = groups[["group"]].copy()

    # strict alignment: keep only overlapping samples
    df = exp.join(grp, how="inner")
    if df.empty:
        raise ValueError("No overlap between exposures.index and groups.index.")

    # sanity: no missing group values after join
    if df["group"].isna().any():
        missing = df.index[df["group"].isna()].tolist()
        raise ValueError(
            f"Missing group labels for samples: {missing[:10]}{'...' if len(missing) > 10 else ''}"
        )

    return df

def plot_one_panel(
    df_panel: pd.DataFrame,
    title: str,
    sig_cols: list[str],
    ylabel: str,
    force_ylim_01: bool = False,
    stacked: bool = True,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    group_labels = df_panel["group"].to_numpy()
    exp_panel = df_panel[sig_cols]

    n_samples = exp_panel.shape[0]
    if figsize is None:
        fig_w = max(6.0, 0.4 * n_samples)
        fig_h = 4.0
        _figsize = (fig_w, fig_h)
    else:
        _figsize = figsize

    fig, ax = plt.subplots(figsize=_figsize)

    # x positions with gaps between groups
    if n_samples == 0:
        x = np.arange(0, dtype=float)
    else:
        x = np.zeros(n_samples, dtype=float)
        x[0] = 0.0
        gap = 1.0
        for i in range(1, n_samples):
            x[i] = x[i - 1] + 1.0
            if group_labels[i] != group_labels[i - 1]:
                x[i] += gap

    if stacked:
        bottom = np.zeros(n_samples)
        for col in sig_cols:  # consistent signature order everywhere
            vals = exp_panel[col].to_numpy()
            ax.bar(x, vals, bottom=bottom, label=str(col))
            bottom += vals
    else:
        n_sig = len(sig_cols)
        width = 0.8 / max(n_sig, 1)
        for i, col in enumerate(sig_cols):
            vals = exp_panel[col].to_numpy()
            ax.bar(x + i * width, vals, width=width, label=str(col))
        if n_samples > 0:
            ax.set_xlim(x.min() - 0.5, x.max() + 0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(exp_panel.index, rotation=90, fontsize="x-small")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize="small", title="Signatures")

    if force_ylim_01:
        ax.set_ylim(0, 1)

    # group labels: start text at the group midpoint (or at first bar if only one group)
    if n_samples > 0:
        starts = [0] + [i for i in range(1, n_samples) if group_labels[i] != group_labels[i - 1]]
        ends = [s - 1 for s in starts[1:]] + [n_samples - 1]
        n_groups = len(starts)

        for s, e in zip(starts, ends):
            x_start = x[s] if n_groups == 1 else (x[s] + x[e]) / 2.0

            ax.text(
                x_start,
                1.06,  # a bit higher than 1.02 to avoid title/legend collisions
                str(group_labels[s]),
                ha="left",      # start text at x_start
                va="bottom",
                rotation=30,
                fontsize="x-small",
                transform=ax.get_xaxis_transform(),  # x in data, y in axes coords
                clip_on=False,
            )

        # Reserve a bit of room at the top for the group labels
        fig.subplots_adjust(top=0.85)

    fig.tight_layout()
    return fig



def signature_sort_key(name: str) -> int:
    m = re.match(r"^Signature_(\d+)$", str(name))
    return int(m.group(1)) if m else 10**9  # non-matching go last


def plot_exposures(
    result: NMFResult,
    *,
    stacked: bool = True,
    figsize: tuple[float, float] | None = None,
    max_samples_per_fig: int | None = None,
    plot: Literal["both", "absolute", "proportion"] = "both",
) -> dict[str, list[plt.Figure]]:
    """
    Plot sample exposures from an :class:`NMFResult`.

    This function generates exposure visualizations while enforcing:

    - Sample order is identical across all panels/plots (absolute, proportion,
      and per-signature views use the same ordering rule).
    - Signature stacking/order is consistent (uses ``result.exposures.columns``).

    Expected inputs:

    - ``result.exposures``: :class:`pandas.DataFrame` with samples as index and
      signatures as columns.
    - ``result.groups``: :class:`pandas.DataFrame` with a ``"group"`` column
      (optional), used to determine sample ordering and/or grouping.

    Parameters
    ----------
    result : NMFResult
        Output of :func:`run_nmf` containing exposures and optional grouping
        information.
    stacked : bool, optional
        If ``True``, plot stacked bar charts. If ``False``, plot grouped (side-by-side)
        bars where applicable. Default is ``True``.
    figsize : tuple[float, float] or None, optional
        Figure size passed to matplotlib. If ``None``, a default size is chosen
        based on the number of samples and signatures.
    max_samples_per_fig : int or None, optional
        Maximum number of samples to show per figure. If provided, samples are
        split across multiple figures. If ``None``, all samples are plotted in a
        single figure (may be large).
    plot : {"both", "absolute", "proportion"}, optional
        Which exposure views to generate:

        - ``"absolute"``: plot raw exposure values.
        - ``"proportion"``: plot exposures normalized to sum to 1 per sample.
        - ``"both"``: generate both absolute and proportional plots.

        Default is ``"both"``.

    Returns
    -------
    dict[str, list[matplotlib.figure.Figure]]
        Dictionary mapping plot type to a list of created figures. Keys depend on
        ``plot`` (e.g. ``"absolute"``, ``"proportion"``).
    """
    exp = result.exposures
    if not isinstance(exp, pd.DataFrame):
        raise TypeError("result.exposures must be a pandas.DataFrame.")

    # ---- build one table: exposures + group ----
    groups = getattr(result, "groups", None)
    if isinstance(groups, pd.DataFrame) and (not groups.empty) and ("group" in groups.columns):
        df = exp.join(groups[["group"]], how="inner")
        if df.empty:
            raise ValueError("No overlap between exposures.index and groups.index.")
    else:
        df = exp.copy()
        df["group"] = "1"

    # ---- one deterministic sample order for ALL plots ----
    sig_cols = [c for c in df.columns if c != "group"]
    sig_cols = sorted(sig_cols, key=signature_sort_key, reverse=False)
    df["_total"] = df[sig_cols].sum(axis=1)
    df = df.sort_values(by=["group", "_total"], ascending=[True, False], kind="mergesort").drop(
        columns=["_total"]
    )

    # ---- proportional exposures (same row order) ----
    df_prop = df.copy()
    totals = df_prop[sig_cols].sum(axis=1).replace(0, np.nan)
    df_prop[sig_cols] = df_prop[sig_cols].div(totals, axis=0).fillna(0.0)
    df_prop = df_prop.sort_values(
        by=["group"] + sig_cols, ascending=[True] + [False] * len(sig_cols), kind="mergesort"
    )
    # ---- chunking (multiple figures if many samples) ----
    n = df.shape[0]
    if max_samples_per_fig is None or n <= max_samples_per_fig:
        chunks = [np.arange(n)]
    else:
        chunks = [
            np.arange(s, min(n, s + max_samples_per_fig)) for s in range(0, n, max_samples_per_fig)
        ]

    figs: dict[str, list[plt.Figure]] = {}

    if plot in ("both", "absolute"):
        figs["absolute"] = [
            plot_one_panel(
                df.iloc[idx],
                "Sample exposures (absolute)",
                sig_cols,
                "Exposure",
                stacked=stacked,
                figsize=figsize,
            )
            for idx in chunks
        ]

    if plot in ("both", "proportion"):
        figs["proportion"] = [
            plot_one_panel(
                df_prop.iloc[idx],
                "Sample exposures (proportions)",
                sig_cols,
                "Proportion",
                force_ylim_01=True,
                stacked=stacked,
                figsize=figsize,
            )
            for idx in chunks
        ]

    return figs
