import hashlib
import json
import os
import re
import zipfile
from pathlib import Path

import pandas as pd
import pytest

from str_mut_signatures import parse_vcf_files, build_mutation_matrix


def file_hash(path: str) -> str:
    """Calculate MD5 hash of a file."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


@pytest.fixture(scope="session")
def vcf_dir(data_dir, tmp_path_factory):
    """
    Unpack tests/data/vcfs/test_input.zip into a temporary directory
    and return the directory containing VCF files.
    """
    zip_path = os.path.join(data_dir, "vcfs", "test_input.zip")
    assert os.path.isfile(zip_path), f"Missing test input zip: {zip_path}"

    tmp_dir = tmp_path_factory.mktemp("vcfs")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(tmp_dir)

    inner = list(Path(tmp_dir).iterdir())
    if len(inner) == 1 and inner[0].is_dir():
        return str(inner[0])
    return str(tmp_dir)


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
        dict(ru="length", ref_length=True, change=True),
        r"^LEN\d+_\d+_[+-]\d+$",
        "bf7b088df766afae97ee38bd3ec59193",
    ),
    (
        "matrix_ru_seq_ref_change",
        dict(ru="ru", ref_length=True, change=True),
        r"^[^_]+_\d+_[+-]\d+$",
        "adddbd6d20b305dd7c1bcdfa73179eb1",
    ),
    (
        "matrix_no_ru_ref_change",
        dict(ru=None, ref_length=True, change=True),
        r"^\d+_[+-]\d+$",
        "0f1b42c745fae72b61151b66222d29f9",
    ),
    (
        "matrix_ru_length_no_change",
        dict(ru="length", ref_length=True, change=False),
        r"^LEN\d+_\d+$",
        "d14d3ba6ce19e548b6c175f518ac5dac",
    ),
    (
        "matrix_ru_seq_change_only",
        dict(ru="ru", ref_length=False, change=True),
        r"^[^_]+_[+-]\d+$",
        "5e397b7760641ef785e490030ac38bac",
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


class TestBuildMutationMatrix:
    @staticmethod
    def _assert_matrix_basic(matrix: pd.DataFrame, name: str):
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

        update_hashes = request.config.getoption("--update-vcf-hashes")
        regex = re.compile(pattern)

        # Build matrix
        matrix = build_mutation_matrix(mutations_df, **kwargs)
        self._assert_matrix_basic(matrix, name)

        # Column name pattern check
        for col in matrix.columns:
            assert regex.match(col), f"{name}: unexpected column name '{col}'"

        # Save matrix
        os.makedirs(output_dir, exist_ok=True)
        out_path = os.path.join(output_dir, f"{name}.tsv")
        matrix.to_csv(out_path, sep="\t")

        actual_hash = file_hash(out_path)

        if update_hashes:
            # Print so you can copy into MATRIX_CASES
            print(f"{name}: {actual_hash}")
        else:
            assert actual_hash == expected_hash, (
                f"{name} hash mismatch:\n"
                f"  expected: {expected_hash}\n"
                f"  actual:   {actual_hash}"
            )
