from __future__ import annotations

import hashlib
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from str_mut_signatures import (
    build_mutation_matrix,
    filter_mutation_matrix,
    load_nmf_result,
    parse_vcf_files,
    plot_exposures,
    plot_pca_samples,
    plot_signatures,
    process_vcf_to_rows,
    project_onto_signatures,
    run_nmf,
    save_nmf_result,
)

# ----------------------------------------------------------------------
# Helpers: hash manifest of everything under a root directory
# ----------------------------------------------------------------------


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

        # Skip dynamic metadata with timestamps
        if path.name == "metadata.json":
            continue

        rel = path.relative_to(root).as_posix()
        h = sha256_file(path)
        entries.append((rel, h))

    lines = [f"{rel}\t{h}" for rel, h in entries]
    return "\n".join(lines) + ("\n" if lines else "")


# ----------------------------------------------------------------------
# Core pipeline helper
# ----------------------------------------------------------------------


def to_figures(obj):
    """
    Normalize various possible return types from plotting functions to
    a flat list of matplotlib Figures.

    Accepts:
      - Figure
      - Axes (uses .figure)
      - sequence of Figures/Axes
    """
    if isinstance(obj, Figure):
        return [obj]
    if isinstance(obj, Axes):
        return [obj.figure]
    if isinstance(obj, Sequence):
        figs = []
        for x in obj:
            if isinstance(x, Figure):
                figs.append(x)
            elif isinstance(x, Axes):
                figs.append(x.figure)
        return figs
    # fallback: anything with savefig()
    if hasattr(obj, "savefig"):
        return [obj]
    return []


def run_full_pipeline(vcf_dir: str, output_dir: str) -> Path:
    """
    Run the full end-to-end pipeline ...
    (docstring unchanged)
    """
    pipeline_dir = Path(output_dir) / "integration_pipeline"
    pipeline_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1) Parse all VCF files into mutations dataframe
    # ------------------------------------------------------------------
    mutations = parse_vcf_files(vcf_dir)
    assert isinstance(mutations, pd.DataFrame)
    assert not mutations.empty

    (pipeline_dir / "mutations.tsv").parent.mkdir(parents=True, exist_ok=True)
    mutations.to_csv(pipeline_dir / "mutations.tsv", sep="\t")

    # ------------------------------------------------------------------
    # 2) Build mutation matrix (samples x features)
    # ------------------------------------------------------------------
    matrix = build_mutation_matrix(
        mutations,
        ru_length=True,
        ru=None,
        ref_length=True,
        change=True,
    )
    assert isinstance(matrix, pd.DataFrame)
    assert not matrix.empty

    matrix.to_csv(pipeline_dir / "matrix_raw.tsv", sep="\t")

    # ------------------------------------------------------------------
    # 3) Filter matrix with all filter methods
    # ------------------------------------------------------------------
    matrix_manual, summary_manual = filter_mutation_matrix(
        matrix,
        feature_method="manual",
        min_feature_total=1,
        min_samples_with_feature=1,
        min_sample_total=1,
    )
    assert not matrix_manual.empty
    assert summary_manual.feature_threshold_used == 1
    matrix_manual.to_csv(pipeline_dir / "matrix_manual.tsv", sep="\t")

    matrix_elbow, summary_elbow = filter_mutation_matrix(
        matrix,
        feature_method="elbow",
        min_samples_with_feature=1,
        min_sample_total=1,
    )
    assert isinstance(matrix_elbow, pd.DataFrame)
    matrix_elbow.to_csv(pipeline_dir / "matrix_elbow.tsv", sep="\t")

    matrix_pct, summary_pct = filter_mutation_matrix(
        matrix,
        feature_method="percentile",
        feature_percentile=0.5,
        min_samples_with_feature=1,
        min_sample_total=1,
    )
    assert isinstance(matrix_pct, pd.DataFrame)
    matrix_pct.to_csv(pipeline_dir / "matrix_percentile.tsv", sep="\t")

    nmf_input = matrix_manual
    assert nmf_input.shape[0] >= 2
    assert nmf_input.shape[1] >= 2

    # ------------------------------------------------------------------
    # 4) Run NMF on filtered matrix
    # ------------------------------------------------------------------
    nmf_res = run_nmf(nmf_input, n_signatures=2, random_state=0)

    assert nmf_res.signatures.shape[1] == 2
    assert nmf_res.exposures.shape[1] == 2
    assert list(nmf_res.signatures.index) == list(nmf_input.columns)
    assert list(nmf_res.exposures.index) == list(nmf_input.index)

    nmf_res.signatures.to_csv(pipeline_dir / "nmf_signatures.tsv", sep="\t")
    nmf_res.exposures.to_csv(pipeline_dir / "nmf_exposures.tsv", sep="\t")

    # ------------------------------------------------------------------
    # 5) Plot signatures, exposures, PCA  (UPDATED)
    # ------------------------------------------------------------------
    sig_obj = plot_signatures(nmf_res, top_n=10)
    exp_obj = plot_exposures(nmf_res, max_samples_per_fig=200)
    coords, var_ratio, pca_obj = plot_pca_samples(nmf_res, title="PCA of exposures")

    # Basic PCA checks unchanged
    assert coords.shape[1] == 2
    assert len(var_ratio) == 2
    coords.to_csv(pipeline_dir / "pca_coords.tsv", sep="\t")
    pd.Series(var_ratio, index=[f"PC{i + 1}" for i in range(len(var_ratio))]).to_csv(
        pipeline_dir / "pca_explained_variance.tsv",
        sep="\t",
    )

    # Normalize to lists of Figures
    sig_figs = to_figures(sig_obj)
    exp_figs = to_figures(exp_obj["absolute"])
    pca_figs = to_figures(pca_obj)

    # Save only the first figure of each type for the integration test
    if sig_figs:
        sig_figs[0].savefig(pipeline_dir / "plot_signatures.png", dpi=150)
    if not exp_figs:
        raise AssertionError("plot_exposures did not return any figures")

    # Save all exposure figures
    if len(exp_figs) == 1:
        # keep old filename for backward compatibility
        exp_figs[0].savefig(pipeline_dir / "plot_exposures.png", dpi=150)
    else:
        for i, f in enumerate(exp_figs, start=1):
            f.savefig(pipeline_dir / f"plot_exposures_{i}.png", dpi=150)

    if pca_figs:
        pca_figs[0].savefig(pipeline_dir / "plot_pca.png", dpi=150)

    # Close all figures to avoid leaking resources
    for fig in sig_figs + exp_figs + pca_figs:
        plt.close(fig)

    # ------------------------------------------------------------------
    # 6) Save & reload NMF result (standard format)
    # ------------------------------------------------------------------
    outdir = Path(output_dir) / "integration_nmf"
    save_nmf_result(nmf_res, outdir)

    assert (outdir / "signatures.tsv").is_file()
    assert (outdir / "exposures.tsv").is_file()
    assert (outdir / "metadata.json").is_file()

    nmf_loaded = load_nmf_result(outdir)

    assert nmf_loaded.signatures.shape == nmf_res.signatures.shape
    assert nmf_loaded.exposures.shape == nmf_res.exposures.shape
    assert list(nmf_loaded.signatures.index) == list(nmf_res.signatures.index)
    assert list(nmf_loaded.exposures.index) == list(nmf_res.exposures.index)

    # ------------------------------------------------------------------
    # 7) Project a *new* VCF onto learned signatures
    # ------------------------------------------------------------------
    vcf_files = sorted(Path(vcf_dir).glob("*.vcf*"))
    assert vcf_files, f"No VCF files found in {vcf_dir}"
    first_vcf = vcf_files[0]

    rows = process_vcf_to_rows(str(first_vcf))
    new_mut_df = pd.DataFrame(rows)
    assert isinstance(new_mut_df, pd.DataFrame)
    assert not new_mut_df.empty

    new_mut_df.to_csv(pipeline_dir / "new_mutations_single_vcf.tsv", sep="\t")

    new_matrix = build_mutation_matrix(
        new_mut_df,
        ru_length=True,
        ru=None,
        ref_length=True,
        change=True,
    )
    assert isinstance(new_matrix, pd.DataFrame)
    assert not new_matrix.empty

    new_matrix.to_csv(pipeline_dir / "new_matrix_single_vcf.tsv", sep="\t")

    exposures_new = project_onto_signatures(
        new_matrix=new_matrix,
        signatures=nmf_loaded.signatures,
        method="nnls",
    )

    assert list(exposures_new.columns) == list(nmf_loaded.signatures.columns)
    assert exposures_new.shape[0] == new_matrix.shape[0]
    assert np.all(exposures_new.to_numpy() >= -1e-10)

    exposures_new.to_csv(pipeline_dir / "new_exposures_single_vcf.tsv", sep="\t")

    return pipeline_dir


class TestFullPipelineIntegration:
    @pytest.mark.integration
    def test_full_pipeline_core(self, vcf_dir: str, output_dir: str):
        """
        Full functional integration test that must pass on all Python versions.
        Runs the entire pipeline and checks core invariants, but does NOT
        enforce exact file hashes.
        """
        pipeline_dir = run_full_pipeline(vcf_dir, output_dir)
        # Sanity: some key files exist
        assert (pipeline_dir / "mutations.tsv").is_file()
        assert (pipeline_dir / "matrix_raw.tsv").is_file()
        assert (pipeline_dir / "nmf_signatures.tsv").is_file()
        assert (pipeline_dir / "nmf_exposures.tsv").is_file()

    # @pytest.mark.integration
    # @pytest.mark.skipif(
    #     sys.version_info < (3, 9),
    #     reason="NMF numerical differences on Python <3.9 change snapshot hashes",
    # )
    # def test_full_pipeline_snapshot(self, vcf_dir: str, output_dir: str, data_dir: str):
    #     """
    #     Snapshot/hash test: compare the manifest of all files under
    #     the pipeline directory to a stored golden file.

    #     This is only enforced on Python >= 3.9, where NMF numerics
    #     are stable enough for exact hash comparison.
    #     """
    #     pipeline_dir = run_full_pipeline(vcf_dir, output_dir)
    #     manifest = build_hash_manifest(pipeline_dir)

    #     gold_path = (
    #         Path(data_dir) / "test_full_pipeline_from_vcf_to_projection_and_snapshot.txt"
    #     )

    #     # If it's somehow missing, bootstrap it once
    #     if not gold_path.exists():
    #         gold_path.write_text(manifest)
    #         pytest.skip(f"Bootstrapped golden manifest at {gold_path}")

    #     expected = gold_path.read_text()
    #     assert manifest == expected
