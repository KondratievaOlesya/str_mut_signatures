from __future__ import annotations

import hashlib
import os
import re

import pandas as pd
import pytest

from str_mut_signatures.extract_tally.extract_mutations import parse_vcf_files
from str_mut_signatures.extract_tally.tally import (
    build_mutation_matrix,
    motif_base_class,
    validate_mutations_data,
)


def file_hash(path: str) -> str:
    """Calculate MD5 hash of a file."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


@pytest.fixture(scope="session")
def mutations_df(vcf_dir):
    """
    Run parse_vcf_files on the prepared test VCFs.
    """
    df = parse_vcf_files(vcf_dir, n_jobs=1)

    assert isinstance(df, pd.DataFrame)
    assert not df.empty, "Parsed mutations DataFrame is empty; check test_input.zip content."

    expected_cols = {
        "sample",
        "tmp_id",
        "tumor_allele_a",
        "tumor_allele_b",
        "normal_allele_a",
        "normal_allele_b",
        "end",
        "period",
        "ref",
        "motif",
    }
    missing = expected_cols - set(df.columns)
    assert not missing, f"Mutations DataFrame missing required columns: {missing}"

    return df


# Matrix configs: name, kwargs for build_mutation_matrix, column name regex, expected hash
MATRIX_CASES = [
    (
        "matrix_ru_length_ref_change",
        {"ru_length": True, "ru": None, "ref_length": True, "change": True},
        r"^LEN\d+_\d+_[+-]\d+$",
        "3ca8aaf98fe1be026aa6c6e57c68fced",
    ),
    (
        "matrix_ru_seq_ref_change",
        {"ru_length": False, "ru": "ru", "ref_length": True, "change": True},
        r"^[^_]+_\d+_[+-]\d+$",
        "cabc9683eaa635f5487e04c056850ee8",
    ),
    (
        "matrix_no_ru_ref_change",
        {"ru_length": False, "ru": None, "ref_length": True, "change": True},
        r"^\d+_[+-]\d+$",
        "9150ab9572cea9faf2939559f37cc254",
    ),
    (
        "matrix_ru_length_no_change",
        {"ru_length": True, "ru": None, "ref_length": True, "change": False},
        r"^LEN\d+_\d+$",
        "3cdb22ae5fbeb6a9ec7b654a72326842",
    ),
    (
        "matrix_ru_seq_change_only",
        {"ru_length": False, "ru": "ru", "ref_length": False, "change": True},
        r"^[^_]+_[+-]\d+$",
        "4579bf828bc9de923e815bf673e70cb1",
    ),
]


@pytest.fixture(params=MATRIX_CASES, scope="session")
def matrix_case(request, mutations_df, output_dir):
    """
    Provides (name, kwargs, pattern, expected_hash, mutations_df, output_dir, request)
    for each matrix configuration.
    """
    name, kwargs, pattern, expected_hash = request.param
    return name, kwargs, pattern, expected_hash, mutations_df, output_dir, request


class TestValidateMutationsData:
    def test_valid_with_motif_column(self):
        df = pd.DataFrame(
            {
                "sample": ["s1"],
                "normal_allele_a": ["10"],
                "normal_allele_b": ["10"],
                "tumor_allele_a": ["10"],
                "tumor_allele_b": ["11"],
                "motif": ["A"],
            }
        )

        motif_col = validate_mutations_data(df)
        assert motif_col == "motif"

    def test_valid_with_RU_column_and_genotype_separator(self):
        df = pd.DataFrame(
            {
                "sample": ["s1"],
                "normal_allele_a": ["10"],
                "normal_allele_b": ["10"],
                "tumor_allele_a": ["10"],
                "tumor_allele_b": ["11"],
                "RU": ["AT"],
                "genotype_separator": ["|"],
            }
        )

        motif_col = validate_mutations_data(df)
        assert motif_col == "RU"

    def test_missing_required_columns_raises(self):
        df = pd.DataFrame({"sample": ["s1"], "motif": ["A"]})
        with pytest.raises(ValueError) as excinfo:
            validate_mutations_data(df)
        msg = str(excinfo.value)
        assert "missing required columns" in msg

    def test_missing_motif_and_RU_raises(self):
        df = pd.DataFrame(
            {
                "sample": ["s1"],
                "normal_allele_a": ["10"],
                "normal_allele_b": ["10"],
                "tumor_allele_a": ["10"],
                "tumor_allele_b": ["11"],
            }
        )
        with pytest.raises(ValueError) as excinfo:
            validate_mutations_data(df)
        msg = str(excinfo.value)
        assert "must contain 'motif' or 'RU'" in msg

    def test_non_dataframe_raises_type_error(self):
        with pytest.raises(TypeError):
            validate_mutations_data("not a dataframe")  # type: ignore[arg-type]

class TestMotifBaseClass:
    def test_at_only_simple(self):
        assert motif_base_class("A") == "AT_only"
        assert motif_base_class("T") == "AT_only"
        assert motif_base_class("AT") == "AT_only"
        assert motif_base_class("tTaA") == "AT_only"

    def test_gc_only_simple(self):
        assert motif_base_class("G") == "GC_only"
        assert motif_base_class("C") == "GC_only"
        assert motif_base_class("GC") == "GC_only"
        assert motif_base_class("cCgG") == "GC_only"

    def test_mixed_with_at_and_gc(self):
        assert motif_base_class("AC") == "mixed"
        assert motif_base_class("AGT") == "mixed"
        assert motif_base_class("CGT") == "mixed"
        assert motif_base_class("ATGC") == "mixed"

    def test_missing_or_empty_motif_returns_na(self):
        assert pd.isna(motif_base_class(None))
        assert pd.isna(motif_base_class(pd.NA))
        assert pd.isna(motif_base_class(""))

    def test_invalid_characters_return_na(self):
        assert pd.isna(motif_base_class("N"))
        assert pd.isna(motif_base_class("ATN"))
        assert pd.isna(motif_base_class("AT-"))
        assert pd.isna(motif_base_class("123"))


class TestBuildMutationMatrix:
    def test_phased_matrix_two_alleles(self):
        """
        Phased genotypes are compared allele-by-allele.

        One locus with:
          normal: 10,10
          tumor : 10,11

        -> allele A delta = 0 (filtered out)
        -> allele B delta = +1
        -> one feature: LEN1_10_+1 with count 1
        """
        df = pd.DataFrame(
            {
                "sample": ["s1"],
                "normal_allele_a": ["10"],
                "normal_allele_b": ["10"],
                "tumor_allele_a": ["10"],
                "tumor_allele_b": ["11"],
                "motif": ["A"],
                "genotype_separator": ["|"],
            }
        )

        mat = build_mutation_matrix(df, ru_length=True, ru=None, ref_length=True, change=True)

        assert mat.shape == (1, 1)
        assert list(mat.index) == ["s1"]
        assert list(mat.columns) == ["LEN1_10_+1"]
        assert mat.loc["s1", "LEN1_10_+1"] == 1

    def test_unphased_matrix_uses_sorted_pairing(self):
        """
        Unphased genotypes are compared after sorting allele pairs, not by summing.

        One locus with:
          normal: 10,10
          tumor : 10,11

        sorted normal = (10,10)
        sorted tumor  = (10,11)

        -> allele A delta = 0 (filtered out)
        -> allele B delta = +1
        -> one feature: LEN1_10_+1 with count 1
        """
        df = pd.DataFrame(
            {
                "sample": ["s1"],
                "normal_allele_a": ["10"],
                "normal_allele_b": ["10"],
                "tumor_allele_a": ["10"],
                "tumor_allele_b": ["11"],
                "motif": ["A"],
                "genotype_separator": ["/"],
            }
        )

        mat = build_mutation_matrix(df, ru_length=True, ru=None, ref_length=True, change=True)

        assert mat.shape == (1, 1)
        assert list(mat.index) == ["s1"]
        assert list(mat.columns) == ["LEN1_10_+1"]
        assert mat.loc["s1", "LEN1_10_+1"] == 1

    def test_unphased_matrix_does_not_collapse_to_total_sum(self):
        """
        Unphased comparison should preserve genotype changes that would be lost
        if allele sums were used.

        normal: 10,12
        tumor : 11,11

        Correct sorted comparison:
          sorted normal = (10,12)
          sorted tumor  = (11,11)
          deltas = +1, -1

        -> two mutation features should be counted
        """
        df = pd.DataFrame(
            {
                "sample": ["s1"],
                "normal_allele_a": ["10"],
                "normal_allele_b": ["12"],
                "tumor_allele_a": ["11"],
                "tumor_allele_b": ["11"],
                "motif": ["A"],
                "genotype_separator": ["/"],
            }
        )

        mat = build_mutation_matrix(df, ru_length=True, ru=None, ref_length=True, change=True)

        assert list(mat.index) == ["s1"]
        assert set(mat.columns) == {"LEN1_10_+1", "LEN1_12_-1"}
        assert mat.loc["s1", "LEN1_10_+1"] == 1
        assert mat.loc["s1", "LEN1_12_-1"] == 1

    def test_at_mode_in_matrix(self):
        """
        RU class mode should classify motifs by base composition.

        Two phased loci for one sample:
        - motif AT -> AT_only
        - motif AC -> mixed

        Each locus has both alleles changed by +1, so each feature count is 2.
        """
        df = pd.DataFrame(
            {
                "sample": ["s1", "s1"],
                "normal_allele_a": ["10", "8"],
                "normal_allele_b": ["10", "8"],
                "tumor_allele_a": ["11", "9"],
                "tumor_allele_b": ["11", "9"],
                "motif": ["AT", "AC"],
                "genotype_separator": ["|", "|"],
            }
        )

        mat = build_mutation_matrix(df, ru_length=False, ru="class", ref_length=True, change=True)

        assert list(mat.index) == ["s1"]
        assert "AT_only_10_+1" in mat.columns
        assert "mixed_8_+1" in mat.columns
        assert mat.loc["s1", "AT_only_10_+1"] == 2
        assert mat.loc["s1", "mixed_8_+1"] == 2

    def test_no_somatic_events_returns_empty_when_features_selected(self):
        """
        When change=True and no allele has non-zero delta, the feature matrix should be empty.
        """
        df = pd.DataFrame(
            {
                "sample": ["s1"],
                "normal_allele_a": ["10"],
                "normal_allele_b": ["10"],
                "tumor_allele_a": ["10"],
                "tumor_allele_b": ["10"],
                "motif": ["A"],
                "genotype_separator": ["|"],
            }
        )

        mat = build_mutation_matrix(df, ru_length=True, ru=None, ref_length=True, change=True)
        assert mat.empty

    def test_non_mutations_filtered_when_no_features_selected(self):
        """
        Even when no feature components are selected, unchanged loci must not be counted.

        Here:
        - first row has no mutation
        - second row has one changed allele

        Expected fallback output:
          mutation_count = 1
        """
        df = pd.DataFrame(
            {
                "sample": ["s1", "s1"],
                "normal_allele_a": ["10", "10"],
                "normal_allele_b": ["10", "10"],
                "tumor_allele_a": ["10", "10"],
                "tumor_allele_b": ["10", "11"],
                "motif": ["A", "A"],
                "genotype_separator": ["|", "|"],
            }
        )

        mat = build_mutation_matrix(
            df,
            ru_length=False,
            ru=None,
            ref_length=False,
            change=False,
        )

        assert list(mat.index) == ["s1"]
        assert list(mat.columns) == ["mutation_count"]
        assert mat.loc["s1", "mutation_count"] == 1

    def test_no_features_selected_counts_only_true_mutations_unphased(self):
        """
        In fallback mutation_count mode, unphased rows should still use sorted pairing.

        normal: 10,12
        tumor : 11,11

        -> sorted pairing gives deltas +1 and -1
        -> mutation_count = 2
        """
        df = pd.DataFrame(
            {
                "sample": ["s1"],
                "normal_allele_a": ["10"],
                "normal_allele_b": ["12"],
                "tumor_allele_a": ["11"],
                "tumor_allele_b": ["11"],
                "motif": ["A"],
                "genotype_separator": ["/"],
            }
        )

        mat = build_mutation_matrix(
            df,
            ru_length=False,
            ru=None,
            ref_length=False,
            change=False,
        )

        assert list(mat.index) == ["s1"]
        assert list(mat.columns) == ["mutation_count"]
        assert mat.loc["s1", "mutation_count"] == 2

    def test_change_false_with_features_still_filters_non_mutations(self):
        """
        When feature components are selected but change=False, unchanged loci should still
        not contribute if the corrected implementation filters non-mutations globally.

        Row 1: unchanged -> should not count
        Row 2: both alleles changed -> should contribute 2 events to LEN1_10
        """
        df = pd.DataFrame(
            {
                "sample": ["s1", "s1"],
                "normal_allele_a": ["10", "10"],
                "normal_allele_b": ["10", "10"],
                "tumor_allele_a": ["10", "11"],
                "tumor_allele_b": ["10", "11"],
                "motif": ["A", "A"],
                "genotype_separator": ["|", "|"],
            }
        )

        mat = build_mutation_matrix(df, ru_length=True, ru=None, ref_length=True, change=False)

        assert list(mat.index) == ["s1"]
        assert list(mat.columns) == ["LEN1_10"]
        assert mat.loc["s1", "LEN1_10"] == 2


class TestBuildMutationMatrixLarge:
    @staticmethod
    def assert_matrix_basic(matrix: pd.DataFrame, name: str):
        assert isinstance(matrix, pd.DataFrame), f"{name}: result is not a DataFrame"
        assert not matrix.empty, f"{name}: matrix is empty"
        assert (matrix.sum(axis=0) > 0).any(), f"{name}: all columns are zero"

    def test_build_and_hash(self, matrix_case):
        (
            name,
            kwargs,
            pattern,
            expected_hash,
            mutations_df,
            output_dir,
            request,
        ) = matrix_case

        regex = re.compile(pattern)

        # Build matrix
        matrix = build_mutation_matrix(mutations_df, **kwargs)
        self.assert_matrix_basic(matrix, name)

        # Column name pattern check
        for col in matrix.columns:
            assert regex.match(col), f"{name}: unexpected column name '{col}'"

        # Save matrix
        os.makedirs(output_dir, exist_ok=True)
        out_path = os.path.join(output_dir, f"{name}.tsv")
        matrix.to_csv(out_path, sep="\t")

        actual_hash = file_hash(out_path)

        assert actual_hash == expected_hash, (
            f"{name} hash mismatch:\n  expected: {expected_hash}\n  actual:   {actual_hash}"
        )
