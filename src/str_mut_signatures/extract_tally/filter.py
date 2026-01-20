from __future__ import annotations

from typing import Literal, NamedTuple

import numpy as np
import pandas as pd


class FilterSummary(NamedTuple):
    """
    Small structured summary of what was filtered.

    Parameters
    ----------
    feature_stats : pandas.DataFrame
        Per-feature filtering statistics.
    sample_stats : pandas.DataFrame
        Per-sample filtering statistics.
    feature_threshold_used : int or None
        Threshold applied for feature filtering, or ``None`` if no feature-level
        threshold was applied.
    sample_threshold_used : int or None
        Threshold applied for sample filtering, or ``None`` if no sample-level
        threshold was applied.

    Attributes
    ----------
    feature_stats : pandas.DataFrame
        Per-feature filtering statistics.
    sample_stats : pandas.DataFrame
        Per-sample filtering statistics.
    feature_threshold_used : int or None
        Feature-level threshold actually used.
    sample_threshold_used : int or None
        Sample-level threshold actually used.
    """

    feature_stats: pd.DataFrame
    sample_stats: pd.DataFrame
    feature_threshold_used: int | None
    sample_threshold_used: int | None


FeatureMethod = Literal["manual", "elbow", "percentile"]


def compute_feature_stats(matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-feature (column) statistics.

    The returned DataFrame is indexed by feature name and contains:

    - ``total_count``: sum over samples
    - ``n_samples_nonzero``: number of samples where feature count > 0
    - ``mean_per_nonzero``: mean count among samples where feature count > 0

    Parameters
    ----------
    matrix : pandas.DataFrame
        Input matrix with samples as rows and features as columns.

    Returns
    -------
    pandas.DataFrame
        Per-feature summary statistics indexed by feature name.
    """
    totals = matrix.sum(axis=0)
    nonzero = (matrix > 0).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_nonzero = totals / nonzero.replace(0, np.nan)

    stats = pd.DataFrame(
        {
            "total_count": totals,
            "n_samples_nonzero": nonzero,
            "mean_per_nonzero": mean_nonzero,
        }
    )
    return stats


def compute_sample_stats(matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-sample (row) statistics.

    The returned DataFrame is indexed by sample identifier and contains:

    - ``total_count``: sum over all features
    - ``n_features_nonzero``: number of features with count > 0

    Parameters
    ----------
    matrix : pandas.DataFrame
        Input matrix with samples as rows and features as columns.

    Returns
    -------
    pandas.DataFrame
        Per-sample summary statistics indexed by sample ID.
    """
    totals = matrix.sum(axis=1)
    nonzero = (matrix > 0).sum(axis=1)

    stats = pd.DataFrame(
        {
            "total_count": totals,
            "n_features_nonzero": nonzero,
        }
    )
    return stats


def elbow_threshold_from_counts(counts: np.ndarray) -> int:
    """
    Estimate a threshold using a simple elbow heuristic.

    Idea:

    - Sort counts descending
    - Consider points ``(i, count_i)`` in 2D
    - Draw a line between the first and last points
    - Find the point with maximum distance from this line (the "elbow")
    - Use its y-value as the threshold

    Parameters
    ----------
    counts : numpy.ndarray
        One-dimensional array of non-negative counts.

    Returns
    -------
    int
        Estimated threshold value (>= 0).
    """
    counts = np.asarray(counts, dtype=float)
    if counts.size == 0:
        return 0

    y = np.sort(counts)[::-1]  # descending
    x = np.arange(len(y), dtype=float)

    # line from first to last
    a = np.array([x[0], y[0]], dtype=float)
    b = np.array([x[-1], y[-1]], dtype=float)
    ab = b - a
    ab_norm = np.linalg.norm(ab)
    if ab_norm == 0:
        # all counts identical
        return int(y[0])

    # distance from each point to the line using 2D cross product magnitude:
    # |a x b| = |a_x * b_y - a_y * b_x|
    distances = []
    for i in range(len(y)):
        p = np.array([x[i], y[i]], dtype=float)
        ap = p - a
        cross_val = ab[0] * ap[1] - ab[1] * ap[0]
        dist = abs(cross_val) / ab_norm
        distances.append(dist)

    idx = int(np.argmax(distances))
    thr = int(y[idx])
    if thr < 0:
        thr = 0
    return thr


def filter_mutation_matrix(
    matrix: pd.DataFrame,
    *,
    # feature-level filtering
    feature_method: FeatureMethod = "manual",
    min_feature_total: int | None = 10,
    min_samples_with_feature: int | None = 3,
    feature_percentile: float = 0.95,
    # sample-level filtering
    min_sample_total: int | None = 0,
) -> tuple[pd.DataFrame, FilterSummary]:
    """
    Filter a mutation count matrix (samples × features) based on simple metrics.

    Parameters
    ----------
    matrix : pandas.DataFrame
        Mutation count matrix with samples as rows and mutation features as columns.

    feature_method : {"manual", "elbow", "percentile"}, optional
        Strategy for choosing a feature-level total-count threshold.

        - ``"manual"``:
            Use ``min_feature_total`` directly.
        - ``"elbow"``:
            Ignore ``min_feature_total`` and use an elbow heuristic based on the
            distribution of feature total counts.
        - ``"percentile"``:
            Ignore ``min_feature_total`` and keep features whose total count is
            >= the given percentile of the distribution.

    min_feature_total : int or None, optional
        Minimal total count across all samples for a feature to be kept
        (only used when ``feature_method="manual"``). If ``None``, no total-count
        threshold is applied.

    min_samples_with_feature : int or None, optional
        Minimal number of samples in which a feature must be non-zero.
        If ``None``, no prevalence threshold is applied.

    feature_percentile : float, optional
        When ``feature_method="percentile"``, features with ``total_count`` >= the
        ``feature_percentile`` quantile of the distribution are kept.
        Must be between 0 and 1.

    min_sample_total : int or None, optional
        Minimal total count per sample to be kept.
        If ``None``, no sample-level filter is applied.

    Returns
    -------
    filtered_matrix : pandas.DataFrame
        Matrix with filtered samples/features.

    summary : FilterSummary
        Structured summary containing:

        - ``feature_stats``: DataFrame of per-feature metrics
        - ``sample_stats``: DataFrame of per-sample metrics
        - ``feature_threshold_used``: int or None
        - ``sample_threshold_used``: int or None

    Raises
    ------
    ValueError
        If ``feature_method`` is not one of ``"manual"``, ``"elbow"``, or
        ``"percentile"``.
    ValueError
        If ``feature_method="percentile"`` and ``feature_percentile`` is not in
        the interval ``[0, 1]``.
    """
    if not isinstance(matrix, pd.DataFrame):
        raise TypeError("matrix must be a pandas.DataFrame")

    # Compute stats up front
    feature_stats = compute_feature_stats(matrix)
    sample_stats = compute_sample_stats(matrix)

    # ----- Feature-level threshold -----
    if feature_method not in ("manual", "elbow", "percentile"):
        raise ValueError("feature_method must be one of 'manual', 'elbow', 'percentile'")

    if feature_method == "manual":
        feature_thr = min_feature_total
    else:
        totals = feature_stats["total_count"].values
        if feature_method == "elbow":
            feature_thr = elbow_threshold_from_counts(totals)
        else:  # "percentile"
            if not (0.0 < feature_percentile <= 1.0):
                raise ValueError("feature_percentile must be in (0, 1].")
            q = float(np.quantile(totals, feature_percentile))
            feature_thr = int(q)

    # Build feature mask
    feature_mask = pd.Series(True, index=matrix.columns)

    if feature_thr is not None:
        feature_mask &= feature_stats["total_count"] >= feature_thr

    if min_samples_with_feature is not None:
        feature_mask &= feature_stats["n_samples_nonzero"] >= min_samples_with_feature

    # ----- Sample-level threshold -----
    sample_thr = min_sample_total
    sample_mask = pd.Series(True, index=matrix.index)
    if sample_thr is not None:
        sample_mask &= sample_stats["total_count"] >= sample_thr

    # ----- Apply masks -----
    filtered = matrix.loc[sample_mask, feature_mask]

    # Recompute stats for filtered matrix (optional but nice)
    feature_stats_filtered = (
        compute_feature_stats(filtered) if not filtered.empty else feature_stats.iloc[0:0]
    )
    sample_stats_filtered = (
        compute_sample_stats(filtered) if not filtered.empty else sample_stats.iloc[0:0]
    )

    summary = FilterSummary(
        feature_stats=feature_stats_filtered,
        sample_stats=sample_stats_filtered,
        feature_threshold_used=feature_thr,
        sample_threshold_used=sample_thr,
    )

    return filtered, summary
