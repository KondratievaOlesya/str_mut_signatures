from __future__ import annotations

import gzip
import textwrap
from pathlib import Path

from str_mut_signatures.extract_tally.validate import (
    VCFValidationResult,
    validate_vcf,
)


class TestValidateVCF:
    # --- helpers ----------------------------------------------------------

    def write_vcf(self, output_dir: str, name: str, content: str) -> Path:
        """
        Write a plain-text .vcf file into the shared output_dir.
        """
        path = Path(output_dir) / name
        path.write_text(textwrap.dedent(content).lstrip("\n"))
        return path

    def write_vcf_gz(self, output_dir: str, name: str, content: str) -> Path:
        """
        Write a gzipped .vcf.gz file into the shared output_dir.
        """
        path = Path(output_dir) / name
        with gzip.open(path, "wt") as fh:
            fh.write(textwrap.dedent(content).lstrip("\n"))
        return path

    # ---------- synthetic VCF tests ----------

    def test_valid_str_annotated_paired_phased(self, output_dir: str):
        """
        VCF with required STR INFO/FORMAT fields, two samples,
        and phased GT using '|'.
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
        ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
        #CHROM POS ID REF ALT QUAL FILTER INFO        FORMAT        NORMAL      TUMOR
        chr1  100 .  A   T   .    PASS RU=A;REF=10  GT:REPCN      0|0:10,10   0|1:10,11
        """
        path = self.write_vcf(output_dir, "phased.vcf", vcf_content)

        res = validate_vcf(path)

        assert isinstance(res, VCFValidationResult)
        assert res.has_str_annotations
        assert res.has_paired_samples
        assert res.normal_sample == "NORMAL"
        assert res.tumor_sample == "TUMOR"
        assert res.genotype_separator == "|"

    def test_valid_str_annotated_paired_unphased(self, output_dir: str):
        """
        VCF with required STR fields and unphased GT using '/'.
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
        ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
        #CHROM POS ID REF ALT QUAL FILTER INFO        FORMAT        NORMAL      TUMOR
        chr1  200 .  A   T   .    PASS RU=A;REF=10  GT:REPCN      0/0:10,10   0/1:10,11
        """
        path = self.write_vcf(output_dir, "unphased.vcf", vcf_content)

        res = validate_vcf(path)

        assert res.has_str_annotations
        assert res.has_paired_samples
        assert res.normal_sample == "NORMAL"
        assert res.tumor_sample == "TUMOR"
        assert res.genotype_separator == "/"

    def test_missing_str_annotations(self, output_dir: str):
        """
        VCF missing RU/REF or REPCN should have has_str_annotations == False.
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=SOMEOTHER,Number=1,Type=String,Description="Not RU/REF">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM POS ID REF ALT QUAL FILTER INFO   FORMAT  NORMAL  TUMOR
        chr1  300 .  A   T   .    PASS .        GT      0/0     0/1
        """
        path = self.write_vcf(output_dir, "missing_str.vcf", vcf_content)

        res = validate_vcf(path)

        assert not res.has_str_annotations
        assert res.has_paired_samples
        assert res.normal_sample == "NORMAL"
        assert res.tumor_sample == "TUMOR"

    def test_single_sample_only(self, output_dir: str):
        """
        VCF with only one sample column should report has_paired_samples == False
        and no normal/tumor names.
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
        ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
        #CHROM POS ID REF ALT QUAL FILTER INFO        FORMAT        ONLYSAMPLE
        chr1  400 .  A   T   .    PASS RU=A;REF=10  GT:REPCN      0|1:10,11
        """
        path = self.write_vcf(output_dir, "single_sample.vcf", vcf_content)

        res = validate_vcf(path)

        assert res.has_str_annotations
        assert not res.has_paired_samples
        assert res.normal_sample is None
        assert res.tumor_sample is None

    def test_missing_gt_format_field(self, output_dir: str):
        """
        VCF without GT in FORMAT (only REPCN) should lead to genotype_separator == None.
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
        ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
        ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
        #CHROM POS ID REF ALT QUAL FILTER INFO        FORMAT   NORMAL    TUMOR
        chr1  500 .  A   T   .    PASS RU=A;REF=10  REPCN    10,10     10,11
        """
        path = self.write_vcf(output_dir, "no_gt.vcf", vcf_content)

        res = validate_vcf(path)

        assert res.has_str_annotations
        assert res.has_paired_samples
        assert res.genotype_separator is None

    def test_gzipped_vcf(self, output_dir: str):
        """
        Ensure .vcf.gz is handled correctly (header + GT detection).
        """
        vcf_content = """
        ##fileformat=VCFv4.3
        ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
        ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
        #CHROM POS ID REF ALT QUAL FILTER INFO        FORMAT        NORMAL      TUMOR
        chr1  600 .  A   T   .    PASS RU=A;REF=10  GT:REPCN      0|0:10,10   0|1:10,11
        """
        path = self.write_vcf_gz(output_dir, "phased.vcf.gz", vcf_content)

        res = validate_vcf(path)

        assert res.has_str_annotations
        assert res.has_paired_samples
        assert res.genotype_separator == "|"

    # ---------- real test files -------------------

    def test_pindel_header_vcf_real_file(self, data_dir: str):
        path = Path(data_dir) / "pindel_header.vcf"
        res = validate_vcf(path)

        assert res.has_str_annotations
        assert res.has_paired_samples
        assert res.genotype_separator == "/"

    def test_mutect2_indel_vcf_no_str_annotations(self, data_dir: str):
        path = Path(data_dir) / "mutec2_indel.vcf.gz"

        res = validate_vcf(path)

        assert not res.has_str_annotations
        assert res.has_paired_samples
        assert res.genotype_separator == "|"

