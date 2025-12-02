from __future__ import annotations

import gzip
import hashlib
import textwrap
from pathlib import Path

import pandas as pd
import pytest

from str_mut_signatures.extract_tally.extract_mutations import (
    parse_info,
    parse_copy_number,
    parse_vcf_files,
    process_vcf_to_rows,
)


def file_hash(path: str) -> str:
    """Calculate MD5 hash of a file."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

class TestHelpers:
    def test_parse_info_basic(self):
        info = parse_info("RU=A;REF=10;SOMETHING=foo")
        assert info["RU"] == "A"
        assert info["REF"] == "10"
        assert info["SOMETHING"] == "foo"

    def test_parse_info_handles_missing_equals(self):
        info = parse_info("FLAG;RU=A")
        # 'FLAG' should be ignored; RU parsed correctly
        assert "FLAG" not in info
        assert info["RU"] == "A"

    def test_parse_copy_number_two_values(self):
        a, b = parse_copy_number("10,11")
        assert a == "10"
        assert b == "11"

    def test_parse_copy_number_single_value(self):
        a, b = parse_copy_number("7")
        assert a == "7"
        assert b == "7"  # homozygous assumption

    def test_parse_copy_number_more_than_two_values(self):
        a, b = parse_copy_number("5,6,7")
        assert a == "."
        assert b == "."


class TestProcessVCFToRows:
    # ---------- helpers ----------

    def write_vcf(self, output_dir: str, name: str, content: str) -> Path:
        path = Path(output_dir) / name
        path.write_text(textwrap.dedent(content).lstrip("\n"))
        return path

    def write_vcf_gz(self, output_dir: str, name: str, content: str) -> Path:
        path = Path(output_dir) / name
        with gzip.open(path, "wt") as fh:
            fh.write(textwrap.dedent(content).lstrip("\n"))
        return path

    # ---------- core behaviour ----------

    def test_process_vcf_basic_phased(self, output_dir: str, capsys):
        """
        Single STR-annotated variant, FILTER=PASS, numeric REPCN.
        Should produce exactly one row with expected fields and phased GT info.
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
        ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
        #CHROM POS ID REF ALT QUAL FILTER INFO                         FORMAT        NORMAL         TUMOR
        chr1  100 .  A   T   .    PASS RU=A;REF=10;PERFECT=TRUE      GT:REPCN      0|0:10,10      0|1:10,11
        """
        path = self.write_vcf(output_dir, "one_variant.vcf", vcf_content)

        rows = process_vcf_to_rows(path, filter_by_pass=True, filter_by_perfect=True)
        assert len(rows) == 1
        row = rows[0]

        # sample name from file base name
        assert row["sample"] == "one_variant"

        assert row["tmp_id"] == "chr1_100"
        assert row["normal_allele_a"] == "10"
        assert row["normal_allele_b"] == "10"
        assert row["tumor_allele_a"] == "10"
        assert row["tumor_allele_b"] == "11"
        assert row["motif"] == "A"
        assert row["ref"] == "10"
        assert row["genotype_separator"] == "|"  # phased

    def test_filter_by_pass_flag(self, output_dir: str):
        """
        When filter_by_pass=True, only FILTER=PASS variants should be kept.
        When filter_by_pass=False, both are kept (if other checks pass).
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
        ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
        #CHROM POS ID REF ALT QUAL FILTER   INFO                         FORMAT        NORMAL        TUMOR
        chr1  100 .  A   T   .    PASS     RU=A;REF=10;PERFECT=TRUE     GT:REPCN      0/0:10,10     0/1:10,11
        chr1  200 .  A   T   .    LowQual  RU=A;REF=12;PERFECT=TRUE     GT:REPCN      0/0:12,12     0/1:12,13
        """
        path = self.write_vcf(output_dir, "filter_pass.vcf", vcf_content)

        # strict: only FILTER=PASS
        rows_pass = process_vcf_to_rows(path, filter_by_pass=True, filter_by_perfect=True)
        assert len(rows_pass) == 1
        assert rows_pass[0]["tmp_id"] == "chr1_100"

        # relaxed: ignore FILTER field
        rows_all = process_vcf_to_rows(path, filter_by_pass=False, filter_by_perfect=True)
        assert len(rows_all) == 2
        tmp_ids = {r["tmp_id"] for r in rows_all}
        assert tmp_ids == {"chr1_100", "chr1_200"}

    def test_filter_by_perfect_flag(self, output_dir: str):
        """
        When filter_by_perfect=True, variants with PERFECT=FALSE should be dropped.
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
        ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
        ##INFO=<ID=PERFECT,Number=1,Type=String,Description="perfect flag">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
        #CHROM POS ID REF ALT QUAL FILTER INFO                              FORMAT        NORMAL        TUMOR
        chr1  100 .  A   T   .    PASS RU=A;REF=10;PERFECT=TRUE            GT:REPCN      0/0:10,10     0/1:10,11
        chr1  200 .  A   T   .    PASS RU=A;REF=12;PERFECT=FALSE           GT:REPCN      0/0:12,12     0/1:12,13
        chr1  300 .  A   T   .    PASS RU=A;REF=14                         GT:REPCN      0/0:14,14     0/1:14,15
        """
        path = self.write_vcf(output_dir, "filter_perfect.vcf", vcf_content)

        # strict PERFECT filter: drop PERFECT=FALSE; line without PERFECT is kept
        rows_strict = process_vcf_to_rows(path, filter_by_pass=True, filter_by_perfect=True)
        tmp_ids_strict = {r["tmp_id"] for r in rows_strict}
        assert tmp_ids_strict == {"chr1_100", "chr1_300"}

        # relaxed PERFECT filter: keep all
        rows_relaxed = process_vcf_to_rows(path, filter_by_pass=True, filter_by_perfect=False)
        tmp_ids_relaxed = {r["tmp_id"] for r in rows_relaxed}
        assert tmp_ids_relaxed == {"chr1_100", "chr1_200", "chr1_300"}

    def test_non_numeric_repcn_is_skipped(self, output_dir: str):
        """
        If any allele REPCN is non-numeric, the variant should be skipped.
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
        ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
        #CHROM POS ID REF ALT QUAL FILTER INFO                         FORMAT        NORMAL          TUMOR
        chr1  100 .  A   T   .    PASS RU=A;REF=10;PERFECT=TRUE      GT:REPCN      0/0:10,10       0/1:10,11
        chr1  200 .  A   T   .    PASS RU=A;REF=12;PERFECT=TRUE      GT:REPCN      0/0:12,NA       0/1:12,13
        """
        path = self.write_vcf(output_dir, "non_numeric.vcf", vcf_content)

        rows = process_vcf_to_rows(path, filter_by_pass=True, filter_by_perfect=True)

        # Only the first variant should be kept
        assert len(rows) == 1
        assert rows[0]["tmp_id"] == "chr1_100"

    def test_invalid_vcf_missing_str_annotations_raises(self, output_dir: str):
        """
        VCF without RU/REF/REPCN should cause process_vcf_to_rows to raise,
        because validate_vcf fails has_str_annotations.
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=SOMEOTHER,Number=1,Type=String,Description="Not RU/REF">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM POS ID REF ALT QUAL FILTER INFO   FORMAT  NORMAL  TUMOR
        chr1  300 .  A   T   .    PASS .        GT      0/0     0/1
        """
        path = self.write_vcf(output_dir, "missing_str.vcf", vcf_content)

        with pytest.raises(ValueError) as excinfo:
            process_vcf_to_rows(path)

        msg = str(excinfo.value)
        assert "missing required STR annotations" in msg

    def test_invalid_vcf_single_sample_raises(self, output_dir: str):
        """
        VCF with only one sample should cause process_vcf_to_rows to raise,
        because validate_vcf has_paired_samples is False.
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
        ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
        #CHROM POS ID REF ALT QUAL FILTER INFO                         FORMAT        ONLYSAMPLE
        chr1  400 .  A   T   .    PASS RU=A;REF=10;PERFECT=TRUE      GT:REPCN      0|1:10,11
        """
        path = self.write_vcf(output_dir, "single_sample.vcf", vcf_content)

        with pytest.raises(ValueError) as excinfo:
            process_vcf_to_rows(path)

        msg = str(excinfo.value)
        assert "must contain at least two samples" in msg


class TestParseVCFFiles:
    def test_parse_vcf_files_on_real_vcf_dir(
        self,
        vcf_dir: str,
        output_dir: str,
    ):
        """
        Integration test on real VCFs:

        - Run parse_vcf_files on vcf_dir.
        - Check basic schema / content.
        - Save to a deterministic TSV.
        - Compare SHA256 hash to a known expected value.
        """
        df = parse_vcf_files(vcf_dir)

        # --- basic content checks ---
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
            "genotype_separator",
        }
        assert expected_cols.issubset(df.columns)

        assert len(df) > 0

        unique_seps = set(df["genotype_separator"].unique())
        assert unique_seps.issubset({"|", "/"})

        # --- write and hash output ---
        out_path = Path(output_dir) / "mutation_information.tsv"
        # be explicit to make hash stable
        df.to_csv(out_path, sep="\t", index=False)

        actual_hash = file_hash(out_path)

        EXPECTED_SHA256 = "4976bd5da202306d2b870154a99ead7f"

        # Temporary assertion so test fails clearly until you plug in hash
        assert (
            actual_hash == EXPECTED_SHA256
        ), f"SHA256 mismatch for {out_path}: {actual_hash} != {EXPECTED_SHA256}"


class TestConSTRainVCFs:
    """
    Integration tests on single real VCFs from conSTRain.
    """
    def test_constrain_vcf_single(
        self,
        data_dir: str,
        tmp_path: Path,
    ):
        """
        This specifically checks that support for FORMAT/REPLEN works
        end-to-end and that variants are parsed correctly.
        """
        vcf_path = Path(data_dir) / "constrain_example.vcf.gz"
        assert vcf_path.exists(), f"Missing test file: {vcf_path}"

        rows = process_vcf_to_rows(
            vcf_path,
            filter_by_pass=False,
            filter_by_perfect=False,
        )

        # Basic checks
        assert len(rows) > 0, "Expected at least one parsed variant from conSTRain VCF"

        df = pd.DataFrame(rows)
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
            "genotype_separator",
        }
        assert expected_cols.issubset(df.columns)

        # All alleles should be numeric strings
        for col in ["tumor_allele_a", "tumor_allele_b", "normal_allele_a", "normal_allele_b"]:
            assert df[col].map(str.isnumeric).all(), f"Non-numeric values in {col} for conSTRain VCF"

        unique_seps = set(df["genotype_separator"].unique())
        assert unique_seps.issubset({"|", "/"})

        # RU/REF must not be empty
        assert (df["motif"] != "").all()
        assert (df["ref"] != "").all()

        # Hash regression for this single file
        out_path = tmp_path / "constrain_mutations.tsv"
        df.to_csv(out_path, sep="\t", index=False)

        actual_hash = file_hash(out_path)
        EXPECTED_MD5 = "0903e4cc10aac5bc58068efc67692ed8"
        assert actual_hash == EXPECTED_MD5
