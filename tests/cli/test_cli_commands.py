from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def build_hash_manifest(root: str | Path) -> str:
    """
    Walk `root` recursively, and build a stable text manifest:

        <relative_path>\t<sha256>

    sorted by path, one per line.

    Skips dynamic files like metadata.json.
    """
    root = Path(root)
    entries: list[tuple[str, str]] = []

    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue

        # ignore dynamic metadata with timestamps
        if path.name == "metadata.json":
            continue

        rel = path.relative_to(root).as_posix()
        h = sha256_file(path)
        entries.append((rel, h))

    lines = [f"{rel}\t{h}" for rel, h in entries]
    return "\n".join(lines) + ("\n" if lines else "")


def file_hash(path: str) -> str:
    """Calculate MD5 hash of a file."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


# Skip all CLI tests if the executable is not available on PATH
pytestmark = pytest.mark.skipif(
    shutil.which("str_mut_signatures") is None,
    reason="str_mut_signatures CLI not found on PATH",
)

# -------------------------------------------------------------------------
# MATRIX CASES FOR `extract`
# -------------------------------------------------------------------------

CLI_MATRIX_CASES = [
    (
        "cli_matrix_ru_length_ref_change",
        # extra CLI args (after --vcf-dir and --out-matrix)
        ["--ru-length", "--ref-length", "--change"],
        r"^LEN\d+_\d+_[+-]\d+$",
        "bf7b088df766afae97ee38bd3ec59193",
    ),
    (
        "cli_matrix_ru_seq_ref_change",
        ["--ru", "ru", "--ref-length", "--change"],
        r"^[^_]+_\d+_[+-]\d+$",
        "adddbd6d20b305dd7c1bcdfa73179eb1",
    ),
    (
        "cli_matrix_no_ru_ref_change",
        ["--ref-length", "--change"],
        r"^\d+_[+-]\d+$",
        "0f1b42c745fae72b61151b66222d29f9",
    ),
    (
        "cli_matrix_ru_length_no_change",
        ["--ru-length", "--ref-length"],
        r"^LEN\d+_\d+$",
        "d14d3ba6ce19e548b6c175f518ac5dac",
    ),
    (
        "cli_matrix_ru_seq_change_only",
        ["--ru", "ru", "--change"],
        r"^[^_]+_[+-]\d+$",
        "5e397b7760641ef785e490030ac38bac",
    ),
]


# -------------------------------------------------------------------------
# Tests: basic CLI behaviour
# -------------------------------------------------------------------------
class TestCLIBasicUsage:
    """Basic CLI usage tests."""

    # @TODO Write checks
    def test_cli_help(self):
        """`str_mut_signatures --help` runs and shows subcommands."""
        result = subprocess.run(["str_mut_signatures", "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        # very loose checks
        assert "extract" in result.stdout
        assert "nmf" in result.stdout

    def test_version_command(self):
        """Test --version flag."""
        result = subprocess.run(["str_mut_signatures", "--version"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "2.0.0" in result.stdout

    def test_no_arguments_fails(self):
        """Test that running without arguments fails."""
        result = subprocess.run(["str_mut_signatures"], capture_output=True, text=True)
        assert result.returncode != 0
        assert "required" in result.stderr.lower() or "error" in result.stderr.lower()


# -------------------------------------------------------------------------
# Tests: `extract` command
# -------------------------------------------------------------------------


class TestCLIExtractCommand:
    """Tests for `extract` command."""

    def test_basic_extract(self, vcf_dir, temp_output_dir):
        """Test basic extract command with motif length + ref length + change."""
        output_matrix = os.path.join(temp_output_dir, "counts_len1.tsv")

        result = subprocess.run(
            [
                "str_mut_signatures",
                "extract",
                "--vcf-dir",
                vcf_dir,
                "--out-matrix",
                output_matrix,
                "--ru-length",
                "--ref-length",
                "--change",
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )

        assert result.returncode == 0, f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"

        # Output file exists and is non-empty
        assert os.path.exists(output_matrix), "Output matrix file was not created"
        assert os.path.getsize(output_matrix) > 0, "Output matrix file is empty"

        # Basic TSV structure checks
        df = pd.read_csv(output_matrix, sep="\t", index_col=0)
        assert isinstance(df, pd.DataFrame)
        assert not df.empty, "Extracted matrix is empty"
        # At least some non-zero columns
        assert (df.sum(axis=0) > 0).any(), "All columns in matrix are zero"
        # Column name pattern: LEN{motif_length}_{ref_length}_{change}
        for col in df.columns:
            assert col.startswith("LEN"), f"Unexpected column name: {col}"
            parts = col.split("_")
            assert len(parts) == 3, f"Unexpected column format: {col}"

    def test_extract_with_full_motif(self, vcf_dir, temp_output_dir):
        """Test extract command with full motif (--ru ru)."""
        output_matrix = os.path.join(temp_output_dir, "counts_ru.tsv")

        result = subprocess.run(
            [
                "str_mut_signatures",
                "extract",
                "--vcf-dir",
                vcf_dir,
                "--out-matrix",
                output_matrix,
                "--ru",
                "ru",
                "--ref-length",
                "--change",
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )

        assert result.returncode == 0, f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"

        assert os.path.exists(output_matrix), "Output matrix file was not created"
        df = pd.read_csv(output_matrix, sep="\t", index_col=0)
        assert not df.empty, "Extracted matrix is empty"

        # Full motif keys: {RU}_{ref_length}_{change}
        for col in df.columns:
            parts = col.split("_")
            assert len(parts) == 3, f"Unexpected column format: {col}"

    def test_missing_vcf_dir_fails(self, temp_output_dir):
        """Test error handling when --vcf-dir points to a non-existent directory."""
        output_matrix = os.path.join(temp_output_dir, "counts.tsv")

        result = subprocess.run(
            [
                "str_mut_signatures",
                "extract",
                "--vcf-dir",
                "/nonexistent/vcfs",
                "--out-matrix",
                output_matrix,
                "--ru-length",
                "--ref-length",
                "--change",
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
        assert (
            "no such file" in result.stderr.lower()
            or "not found" in result.stderr.lower()
            or "error" in result.stderr.lower()
        )

    def test_extract_missing_out_matrix_fails(self, vcf_dir):
        """Test that --vcf-dir without --out-matrix fails."""
        result = subprocess.run(
            [
                "str_mut_signatures",
                "extract",
                "--vcf-dir",
                vcf_dir,
                "--ru-length",
                "--ref-length",
                "--change",
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
        assert "out-matrix" in result.stderr.lower() or "required" in result.stderr.lower()


class TestCLINMFCommand:
    """Tests for `nmf` command."""

    def test_nmf_runs_and_produces_outputs(self, vcf_dir, temp_output_dir):
        """Integration test: extract -> nmf."""
        # 1) Build matrix via extract
        matrix_path = os.path.join(temp_output_dir, "counts_for_nmf.tsv")

        extract_result = subprocess.run(
            [
                "str_mut_signatures",
                "extract",
                "--vcf-dir",
                vcf_dir,
                "--out-matrix",
                matrix_path,
                "--ru-length",
                "--ref-length",
                "--change",
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )

        assert extract_result.returncode == 0, (
            "extract for NMF test failed\n"
            f"STDOUT:\n{extract_result.stdout}\n\n"
            f"STDERR:\n{extract_result.stderr}"
        )
        assert os.path.exists(matrix_path), "Matrix for NMF was not created"

        df = pd.read_csv(matrix_path, sep="\t", index_col=0)
        assert not df.empty, "Matrix for NMF is empty"

        # 2) Run nmf
        nmf_outdir = os.path.join(temp_output_dir, "nmf_results")
        os.makedirs(nmf_outdir, exist_ok=True)
        n_signatures = 3

        nmf_result = subprocess.run(
            [
                "str_mut_signatures",
                "nmf",
                "--matrix",
                matrix_path,
                "--outdir",
                nmf_outdir,
                "--n-signatures",
                str(n_signatures),
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )

        assert nmf_result.returncode == 0, (
            f"nmf CLI failed\nSTDOUT:\n{nmf_result.stdout}\n\nSTDERR:\n{nmf_result.stderr}"
        )

        # Check that some TSV outputs exist
        tsv_files = list(Path(nmf_outdir).glob("*.tsv"))
        assert tsv_files, "No TSV outputs created by NMF command"

        # Sanity: at least one of them should have dimension n_signatures
        has_nsig_dim = False
        for tsv in tsv_files:
            out_df = pd.read_csv(tsv, sep="\t", index_col=0)
            assert not out_df.empty, f"NMF output file {tsv} is empty"
            if out_df.shape[0] == n_signatures or out_df.shape[1] == n_signatures:
                has_nsig_dim = True

        assert has_nsig_dim, (
            f"No NMF TSV output appears to have a dimension equal to n_signatures={n_signatures}"
        )

    def test_nmf_missing_arguments_fail(self, temp_output_dir):
        """Test argument validation for nmf command."""
        dummy_matrix = os.path.join(temp_output_dir, "dummy.tsv")
        Path(dummy_matrix).write_text("sample\tfeat1\ns1\t1\n")

        # Missing --outdir
        result_no_outdir = subprocess.run(
            [
                "str_mut_signatures",
                "nmf",
                "--matrix",
                dummy_matrix,
                "--n-signatures",
                "3",
            ],
            capture_output=True,
            text=True,
        )
        assert result_no_outdir.returncode != 0
        assert (
            "outdir" in result_no_outdir.stderr.lower()
            or "required" in result_no_outdir.stderr.lower()
        )

        # Missing --matrix
        result_no_matrix = subprocess.run(
            [
                "str_mut_signatures",
                "nmf",
                "--outdir",
                temp_output_dir,
                "--n-signatures",
                "3",
            ],
            capture_output=True,
            text=True,
        )
        assert result_no_matrix.returncode != 0
        assert (
            "matrix" in result_no_matrix.stderr.lower()
            or "required" in result_no_matrix.stderr.lower()
        )

    def test_nmf_invalid_matrix_path(self, temp_output_dir):
        """Test error handling for invalid matrix path."""
        nmf_outdir = os.path.join(temp_output_dir, "nmf_results_invalid")
        os.makedirs(nmf_outdir, exist_ok=True)

        result = subprocess.run(
            [
                "str_mut_signatures",
                "nmf",
                "--matrix",
                "/nonexistent/matrix.tsv",
                "--outdir",
                nmf_outdir,
                "--n-signatures",
                "3",
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
        assert "no such file" in result.stderr.lower() or "not found" in result.stderr.lower()


@pytest.mark.integration
class TestCLIIntegration:
    """End-to-end CLI integration tests."""

    def run_cli(self, args: list[str]) -> subprocess.CompletedProcess:
        """
        Helper: run the str_mut_signatures CLI with given args and
        assert successful return code.
        """
        cmd = ["str_mut_signatures"] + args
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, (
            f"CLI failed: {' '.join(cmd)}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
        return result

    # ------------------------------------------------------------------
    # Helper: run full CLI pipeline once and return pipeline_dir
    # ------------------------------------------------------------------
    def run_full_cli_pipeline(self, vcf_dir: str, output_dir: str) -> Path:
        """
        1) extract: VCF dir -> count matrix
        2) filter: matrix -> filtered matrix
        3) nmf    : filtered matrix -> signatures + exposures + metadata.json
        4) project: new matrix -> new exposures

        All outputs are written under:
            pipeline_dir = output_dir / "cli_integration_pipeline"
        """
        pipeline_dir = Path(output_dir) / "cli_integration_pipeline"
        pipeline_dir.mkdir(parents=True, exist_ok=True)

        # Paths for intermediate artifacts
        matrix_raw = pipeline_dir / "matrix_raw.tsv"
        matrix_filtered = pipeline_dir / "matrix_filtered.tsv"
        nmf_dir = pipeline_dir / "nmf_results"
        nmf_dir.mkdir(parents=True, exist_ok=True)
        new_matrix = pipeline_dir / "new_matrix_single_vcf.tsv"
        new_exposures = pipeline_dir / "new_exposures_single_vcf.tsv"

        # ------------------------------------------------------------------
        # 1) extract: from all VCFs
        # ------------------------------------------------------------------
        self.run_cli(
            [
                "extract",
                "--vcf-dir",
                vcf_dir,
                "--out-matrix",
                str(matrix_raw),
                "--ru-length",
                "--ref-length",
                "--change",
            ]
        )
        assert matrix_raw.is_file()

        # ------------------------------------------------------------------
        # 2) filter: feature-level filtering
        # ------------------------------------------------------------------
        self.run_cli(
            [
                "filter",
                "--matrix",
                str(matrix_raw),
                "--out-matrix",
                str(matrix_filtered),
                "--feature-method",
                "manual",
                "--min-feature-total",
                "1",
                "--min-samples-with-feature",
                "1",
                "--min-sample-total",
                "1",
            ]
        )
        assert matrix_filtered.is_file()

        # ------------------------------------------------------------------
        # 3) nmf: run decomposition on filtered matrix
        # ------------------------------------------------------------------
        self.run_cli(
            [
                "nmf",
                "--matrix",
                str(matrix_filtered),
                "--outdir",
                str(nmf_dir),
                "--n-signatures",
                "2",
                "--max-iter",
                "200",
                "--random-state",
                "0",
            ]
        )

        # These are created by save_nmf_result inside the CLI
        sig_path = nmf_dir / "signatures.tsv"
        exp_path = nmf_dir / "exposures.tsv"
        meta_path = nmf_dir / "metadata.json"

        assert sig_path.is_file()
        assert exp_path.is_file()
        assert meta_path.is_file()

        # ------------------------------------------------------------------
        # 4) project: build a "new" matrix and project onto learned signatures
        # ------------------------------------------------------------------
        # Create a new VCF directory with a single VCF to mimic a new cohort
        new_vcf_dir = pipeline_dir / "single_vcf"
        new_vcf_dir.mkdir(parents=True, exist_ok=True)

        vcf_files = sorted(Path(vcf_dir).glob("*.vcf*"))
        assert vcf_files, f"No VCF files found in {vcf_dir}"
        first_vcf = vcf_files[0]

        shutil.copy2(first_vcf, new_vcf_dir / first_vcf.name)

        # Extract for this single-VCF cohort
        self.run_cli(
            [
                "extract",
                "--vcf-dir",
                str(new_vcf_dir),
                "--out-matrix",
                str(new_matrix),
                "--ru-length",
                "--ref-length",
                "--change",
            ]
        )
        assert new_matrix.is_file()

        # Project onto signatures from nmf_dir
        self.run_cli(
            [
                "project",
                "--matrix",
                str(new_matrix),
                "--nmf-dir",
                str(nmf_dir),
                "--out-exposures",
                str(new_exposures),
            ]
        )
        assert new_exposures.is_file()

        return pipeline_dir

    # ------------------------------------------------------------------
    # Tests
    # ------------------------------------------------------------------

    def test_full_pipeline_cli_core(self, vcf_dir: str, output_dir: str):
        """
        Full CLI pipeline that must succeed on all Python versions.
        Does NOT enforce exact file hashes.
        """
        pipeline_dir = self.run_full_cli_pipeline(vcf_dir, output_dir)

        # Basic sanity checks for key files
        assert (pipeline_dir / "matrix_raw.tsv").is_file()
        assert (pipeline_dir / "matrix_filtered.tsv").is_file()
        nmf_dir = pipeline_dir / "nmf_results"
        assert (nmf_dir / "signatures.tsv").is_file()
        assert (nmf_dir / "exposures.tsv").is_file()
        assert (nmf_dir / "metadata.json").is_file()
        assert (pipeline_dir / "new_matrix_single_vcf.tsv").is_file()
        assert (pipeline_dir / "new_exposures_single_vcf.tsv").is_file()

    # @pytest.mark.skipif(
    #     sys.version_info < (3, 9),
    #     reason="NMF numerical differences on Python <3.8 change snapshot hashes",
    # )
    # def test_full_pipeline_cli_snapshot(
    #     self,
    #     vcf_dir: str,
    #     output_dir: str,
    #     data_dir: str,
    # ):
    #     """
    #     Snapshot/hash test for CLI pipeline.
    #     """
    #     pipeline_dir = self.run_full_cli_pipeline(vcf_dir, output_dir)
    #     manifest = build_hash_manifest(pipeline_dir)

    #     gold_path = (
    #         Path(data_dir)
    #         / "test_cli_full_pipeline_from_vcf_to_projection_and_snapshot.txt"
    #     )

    #     if not gold_path.exists():
    #         # one-time bootstrap if you haven't created the golden file yet
    #         gold_path.write_text(manifest)
    #         pytest.skip(f"Bootstrapped CLI golden manifest at {gold_path}")

    #     expected = gold_path.read_text()
    #     assert manifest == expected
