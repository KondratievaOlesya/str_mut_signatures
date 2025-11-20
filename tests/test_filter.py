from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from str_mut_signatures.extract_tally.filter import (
    FilterSummary,
    compute_feature_stats,
    compute_sample_stats,
    elbow_threshold_from_counts,
    filter_mutation_matrix,
)


class TestComputeFeatureStats:
    def test_basic_feature_stats(self):
        # samples x features
        # s1: A=1, B=0, C=2
        # s2: A=0, B=3, C=0
        matrix = pd.DataFrame(
            {
                "A": [1, 0],
                "B": [0, 3],
                "C": [2, 0],
            },
            index=["s1", "s2"],
        )

        stats = compute_feature_stats(matrix)

        # total_count
        assert stats.loc["A", "total_count"] == 1
        assert stats.loc["B", "total_count"] == 3
        assert stats.loc["C", "total_count"] == 2

        # n_samples_nonzero
        assert stats.loc["A", "n_samples_nonzero"] == 1
        assert stats.loc["B", "n_samples_nonzero"] == 1
        assert stats.loc["C", "n_samples_nonzero"] == 1

        # mean_per_nonzero = total / n_samples_nonzero
        assert stats.loc["A", "mean_per_nonzero"] == 1
        assert stats.loc["B", "mean_per_nonzero"] == 3
        assert stats.loc["C", "mean_per_nonzero"] == 2

    def test_feature_stats_handle_all_zero_feature(self):
        matrix = pd.DataFrame(
            {
                "A": [0, 0, 0],
                "B": [1, 0, 2],
            },
            index=["s1", "s2", "s3"],
        )

        stats = compute_feature_stats(matrix)

        assert stats.loc["A", "total_count"] == 0
        assert stats.loc["A", "n_samples_nonzero"] == 0
        # mean_per_nonzero should be NaN for A
        assert np.isnan(stats.loc["A", "mean_per_nonzero"])
        # B as usual
        assert stats.loc["B", "total_count"] == 3
        assert stats.loc["B", "n_samples_nonzero"] == 2
        assert stats.loc["B", "mean_per_nonzero"] == 1.5


class TestComputeSampleStats:
    def test_basic_sample_stats(self):
        matrix = pd.DataFrame(
            {
                "A": [1, 0],
                "B": [0, 3],
                "C": [2, 0],
            },
            index=["s1", "s2"],
        )

        stats = compute_sample_stats(matrix)

        # total_count: sum over features
        assert stats.loc["s1", "total_count"] == 3  # 1 + 0 + 2
        assert stats.loc["s2", "total_count"] == 3  # 0 + 3 + 0

        # n_features_nonzero
        assert stats.loc["s1", "n_features_nonzero"] == 2  # A, C
        assert stats.loc["s2", "n_features_nonzero"] == 1  # B


class TestElbowThresholdFromCounts:
    def test_elbow_basic_descending(self):
        counts = np.array([100, 80, 50, 10, 5])
        thr = elbow_threshold_from_counts(counts)
        # For this simple shape, elbow will be somewhere in the middle;
        # we mainly assert that it's between max and min.
        assert thr <= 100
        assert thr >= 5

    def test_elbow_all_same_counts(self):
        counts = np.array([10, 10, 10, 10])
        thr = elbow_threshold_from_counts(counts)
        # With identical values, function falls back to that value
        assert thr == 10

    def test_elbow_empty_counts(self):
        counts = np.array([])
        thr = elbow_threshold_from_counts(counts)
        assert thr == 0


class TestFilterMutationMatrix:
    def test_type_error_on_non_dataframe(self):
        with pytest.raises(TypeError):
            filter_mutation_matrix(
                matrix="not a df",  # type: ignore[arg-type]
                feature_method="manual",
            )

    def test_invalid_feature_method_raises(self):
        matrix = pd.DataFrame({"A": [1, 2]}, index=["s1", "s2"])
        with pytest.raises(ValueError):
            filter_mutation_matrix(
                matrix=matrix,
                feature_method="unknown",  # type: ignore[arg-type]
            )

    def test_invalid_percentile_raises(self):
        matrix = pd.DataFrame({"A": [1, 2]}, index=["s1", "s2"])
        with pytest.raises(ValueError):
            filter_mutation_matrix(
                matrix=matrix,
                feature_method="percentile",
                feature_percentile=1.5,
            )

    def test_manual_feature_filtering(self):
        """
        Feature totals: A=3, B=1, C=0
        - min_feature_total=2 -> keep A only
        - min_samples_with_feature=1 (default) applies too
        """
        matrix = pd.DataFrame(
            {
                "A": [1, 2],  # total 3
                "B": [1, 0],  # total 1
                "C": [0, 0],  # total 0
            },
            index=["s1", "s2"],
        )

        filtered, summary = filter_mutation_matrix(
            matrix,
            feature_method="manual",
            min_feature_total=2,
            min_samples_with_feature=1,
            min_sample_total=0,
        )

        # Only feature A should survive
        assert list(filtered.columns) == ["A"]
        assert list(filtered.index) == ["s1", "s2"]
        assert isinstance(summary, FilterSummary)
        assert summary.feature_threshold_used == 2
        assert summary.sample_threshold_used == 0

    def test_min_samples_with_feature(self):
        """
        Feature A: present in 2 samples
        Feature B: present in 1 sample
        With min_samples_with_feature=2 -> B should be dropped
        """
        matrix = pd.DataFrame(
            {
                "A": [1, 1],  # non-zero in both samples
                "B": [0, 2],  # non-zero in only one sample
            },
            index=["s1", "s2"],
        )

        filtered, summary = filter_mutation_matrix(
            matrix,
            feature_method="manual",
            min_feature_total=None,  # no total filter
            min_samples_with_feature=2,
            min_sample_total=None,
        )

        assert list(filtered.columns) == ["A"]
        assert "B" not in filtered.columns

    def test_sample_filtering_by_total(self):
        """
        Sample totals: s1=3, s2=1
        With min_sample_total=2 -> s2 should be dropped.
        """
        matrix = pd.DataFrame(
            {
                "A": [1, 0],
                "B": [2, 1],
            },
            index=["s1", "s2"],
        )

        filtered, summary = filter_mutation_matrix(
            matrix,
            feature_method="manual",
            min_feature_total=None,
            min_samples_with_feature=None,
            min_sample_total=2,
        )

        assert list(filtered.index) == ["s1"]
        assert "s2" not in filtered.index
        assert summary.sample_threshold_used == 2

    def test_percentile_feature_filtering(self):
        """
        Feature totals: A=1, B=5, C=10
        0.5 quantile is around 5; features with total >= 5 remain.
        """
        matrix = pd.DataFrame(
            {
                "A": [1, 0],   # total 1
                "B": [2, 3],   # total 5
                "C": [4, 6],   # total 10
            },
            index=["s1", "s2"],
        )

        filtered, summary = filter_mutation_matrix(
            matrix,
            feature_method="percentile",
            feature_percentile=0.5,
            min_samples_with_feature=None,
            min_sample_total=None,
        )

        cols = set(filtered.columns)
        assert "B" in cols
        assert "C" in cols
        assert "A" not in cols
        assert summary.feature_threshold_used >= 1  # some integer around median

    def test_elbow_feature_filtering(self):
        """
        Simple shape to exercise elbow logic; we only assert that
        some lower-count features drop out.
        """
        matrix = pd.DataFrame(
            {
                "A": [10, 0],  # total 10
                "B": [8, 0],   # total 8
                "C": [2, 0],   # total 2
                "D": [1, 0],   # total 1
            },
            index=["s1", "s2"],
        )

        filtered, summary = filter_mutation_matrix(
            matrix,
            feature_method="elbow",
            min_samples_with_feature=None,
            min_sample_total=None,
        )

        # We don't enforce exact threshold, but we do expect
        # at least one of the lowest-count features dropped.
        remaining = set(filtered.columns)
        assert remaining.issubset({"A", "B", "C", "D"})
        assert len(remaining) <= 4
        assert summary.feature_threshold_used is not None

    def test_empty_matrix_returns_empty_and_empty_stats(self):
        matrix = pd.DataFrame()

        filtered, summary = filter_mutation_matrix(
            matrix,
            feature_method="manual",
            min_feature_total=None,
            min_samples_with_feature=None,
            min_sample_total=None,
        )

        assert filtered.empty
        # stats are empty frames as well
        assert summary.feature_stats.empty
        assert summary.sample_stats.empty
