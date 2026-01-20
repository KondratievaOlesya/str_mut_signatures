from __future__ import annotations

import hashlib
import os
import re

import pandas as pd
import pytest

from str_mut_signatures.extract_tally.extract_mutations import parse_vcf_files
from str_mut_signatures.extract_tally.tally import (
    build_mutation_matrix,
    compute_changes_for_row,
    is_phased,
    make_feature,
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
    df = parse_vcf_files(vcf_dir)

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
        "afd366b9975dd59ec6dbc593d7bb8ea8",
    ),
    (
        "matrix_ru_seq_ref_change",
        {"ru_length": False, "ru": "ru", "ref_length": True, "change": True},
        r"^[^_]+_\d+_[+-]\d+$",
        "51c062035ebf6a3c19ac295fa085b104",
    ),
    (
        "matrix_no_ru_ref_change",
        {"ru_length": False, "ru": None, "ref_length": True, "change": True},
        r"^\d+_[+-]\d+$",
        "27b480e60805e209500d3839502dd149",
    ),
    (
        "matrix_ru_length_no_change",
        {"ru_length": True, "ru": None, "ref_length": True, "change": False},
        r"^LEN\d+_\d+$",
        "1d5780166d7a508b847a0a8291079e8c",
    ),
    (
        "matrix_ru_seq_change_only",
        {"ru_length": False, "ru": "ru", "ref_length": False, "change": True},
        r"^[^_]+_[+-]\d+$",
        "10a7c3cf0cd73b68e81744aa746a19d6",
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

        motif_col, has_sep = validate_mutations_data(df)
        assert motif_col == "motif"
        assert not has_sep

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

        motif_col, has_sep = validate_mutations_data(df)
        assert motif_col == "RU"
        assert has_sep

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


class TestIsPhased:
    def test_is_phased_true_for_pipe(self):
        assert is_phased("|")

    def test_is_phased_false_for_slash(self):
        assert not is_phased("/")

    def test_is_phased_false_for_none(self):
        assert not is_phased(None)

    def test_is_phased_false_for_other_string(self):
        assert not is_phased("?")


class TestComputeChangesForRow:
    def test_phased_changes_two_alleles(self):
        row = pd.Series(
            {
                "normal_allele_a": "10",
                "normal_allele_b": "11",
                "tumor_allele_a": "10",
                "tumor_allele_b": "14",
                "genotype_separator": "|",
            }
        )

        out = compute_changes_for_row(row)

        assert out["change_a"] == 0  # 10 - 10
        assert out["change_b"] == 3  # 14 - 11
        assert out["ref_a"] == 10
        assert out["ref_b"] == 11

    def test_unphased_combined_change(self):
        row = pd.Series(
            {
                "normal_allele_a": "10",
                "normal_allele_b": "11",
                "tumor_allele_a": "10",
                "tumor_allele_b": "14",
                "genotype_separator": "/",  # unphased
            }
        )

        out = compute_changes_for_row(row)

        # total_normal = 21, total_tumor = 24, delta = +3
        assert out["change_a"] == 3
        assert pd.isna(out["change_b"])
        assert out["ref_a"] == 21
        assert pd.isna(out["ref_b"])

    def test_no_genotype_separator_treated_as_unphased(self):
        row = pd.Series(
            {
                "normal_allele_a": "5",
                "normal_allele_b": "5",
                "tumor_allele_a": "7",
                "tumor_allele_b": "7",
            }
        )

        out = compute_changes_for_row(row)

        # total_normal = 10, total_tumor = 14, delta = +4
        assert out["change_a"] == 4
        assert pd.isna(out["change_b"])
        assert out["ref_a"] == 10
        assert pd.isna(out["ref_b"])

    def test_non_numeric_alleles_return_na(self):
        row = pd.Series(
            {
                "normal_allele_a": "NA",
                "normal_allele_b": "10",
                "tumor_allele_a": "10",
                "tumor_allele_b": "11",
                "genotype_separator": "|",
            }
        )

        out = compute_changes_for_row(row)

        assert all(pd.isna(out[k]) for k in ["change_a", "change_b", "ref_a", "ref_b"])


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


class TestMakeFeature:
    def test_length_mode_with_ref_and_change(self):
        feat = make_feature(
            motif="A",
            ref=10,
            delta=1,
            ru_length=True,
            ru=None,
            ref_length=True,
            change=True,
        )
        assert feat == "LEN1_10_+1"

    def test_ru_mode_with_ref_and_change(self):
        feat = make_feature(
            motif="AT",
            ref=20,
            delta=-2,
            ru_length=False,
            ru="ru",
            ref_length=True,
            change=True,
        )
        assert feat == "AT_20_-2"

    def test_at_mode_at_rich(self):
        feat = make_feature(
            motif="AT",
            ref=10,
            delta=1,
            ru_length=False,
            ru="class",
            ref_length=True,
            change=True,
        )
        assert feat == "AT_only_10_+1"

    def test_at_mode_non_at_rich(self):
        feat = make_feature(
            motif="AC",
            ref=8,
            delta=2,
            ru_length=False,
            ru="class",
            ref_length=True,
            change=True,
        )
        assert feat == "mixed_8_+2"

    def test_no_change_component_when_change_false(self):
        feat = make_feature(
            motif="A",
            ref=10,
            delta=0,
            ru_length=True,
            ru=None,
            ref_length=True,
            change=False,
        )
        # Only LEN1 and ref length
        assert feat == "LEN1_10"

    def test_delta_zero_dropped_when_change_true(self):
        feat = make_feature(
            motif="A",
            ref=10,
            delta=0,
            ru_length=True,
            ru=None,
            ref_length=True,
            change=True,
        )
        assert pd.isna(feat)

    def test_missing_ref_drops_feature_when_ref_length_true(self):
        feat = make_feature(
            motif="A",
            ref=pd.NA,
            delta=1,
            ru_length=True,
            ru=None,
            ref_length=True,
            change=True,
        )
        assert pd.isna(feat)

    def test_ru_none_ref_false_change_false_results_in_na(self):
        feat = make_feature(
            motif="A",
            ref=10,
            delta=1,
            ru_length=False,
            ru=None,
            ref_length=False,
            change=False,
        )
        # nothing added to parts → NA
        assert pd.isna(feat)

    def test_invalid_ru_raises(self):
        with pytest.raises(ValueError):
            make_feature(
                motif="A",
                ref=10,
                delta=1,
                ru_length=False,
                ru="something",  # type: ignore[arg-type]
                ref_length=True,
                change=True,
            )


class TestBuildMutationMatrix:
    def test_phased_matrix_two_alleles(self):
        """
        Phased: we expect allele-specific events.
        One locus with:
          normal: 10,10
          tumor : 10,11
        → change_a = 0 (dropped), change_b = +1 → 1 event.
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
        cols = list(mat.columns)
        assert cols == ["LEN1_10_+1"]
        assert mat.iloc[0, 0] == 1

    def test_unphased_matrix_combined_event(self):
        """
        Unphased: one combined event per locus.
        One locus with:
          normal: 10,10  (total 20)
          tumor : 10,11  (total 21)
        → combined change +1 at ref 20.
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
        cols = list(mat.columns)
        assert cols == ["LEN1_20_+1"]
        assert mat.iloc[0, 0] == 1

    def test_at_mode_in_matrix(self):
        """
        AT mode: classify motifs into AT_rich / non_AT_rich.
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

        cols = set(mat.columns)
        # First row: AT_only_10_+1 (each allele 10 -> 11 → two events, but both same)
        # Second row: mixed_8_+1 (similar)
        assert "AT_only_10_+1" in cols
        assert "mixed_8_+1" in cols
        # Only one sample "s1"
        assert list(mat.index) == ["s1"]

    def test_no_somatic_events_returns_empty(self):
        """
        When all deltas are zero and change=True, matrix should be empty.
        """
        df = pd.DataFrame(
            {
                "sample": ["s1"],
                "normal_allele_a": ["10"],
                "normal_allele_b": ["10"],
                "tumor_allele_a": ["10"],
                "tumor_allele_b": ["10"],  # no change
                "motif": ["A"],
                "genotype_separator": ["|"],
            }
        )

        mat = build_mutation_matrix(df, ru_length=True, ru=None, ref_length=True, change=True)
        assert mat.empty

    def test_change_false_keeps_all_events(self):
        """
        With change=False, features do not depend on delta and all loci
        that pass numeric checks should contribute.
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

        # With change=False + ref_length=True + ru_length=True:
        # feature is LEN1_10 for each allele/event.
        assert list(mat.index) == ["s1"]
        assert list(mat.columns) == ["LEN1_10"]
        # 2 rows * 2 alleles each = 4 events
        assert mat.iloc[0, 0] == 4


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
