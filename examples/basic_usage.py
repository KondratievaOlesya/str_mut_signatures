"""
basic_usage.py

Example end-to-end workflow for `str_mut_signatures`:

1. Parse STR-annotated paired tumor–normal VCFs from a directory.
2. Build a somatic STR mutation count matrix.
3. Filter the matrix (feature/sample filtering).
4. Run NMF to learn STR mutation signatures.
5. Save and reload the NMF result.
6. Project a "new" VCF onto the learned signatures.
7. Optionally create basic plots (signatures, exposures, PCA).

Adjust paths at the bottom and run:

    python basic_usage.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from str_mut_signatures import (
    build_mutation_matrix,
    compute_pca,
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


def run_full_pipeline(vcf_dir: str | Path, outdir: str | Path, n_signatures: int = 2) -> None:
    """
    End-to-end STR mutation signature analysis.

    Parameters
    ----------
    vcf_dir : str or Path
        Directory with STR-annotated, paired tumor–normal VCF files.
    outdir : str or Path
        Output directory for matrices, NMF results, and plots.
    n_signatures : int
        Number of signatures (rank) for NMF.
    """
    vcf_dir = Path(vcf_dir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1) Parse VCFs -> long-format mutation table
    # ------------------------------------------------------------------
    print(f"[1/7] Parsing VCFs from: {vcf_dir}")
    mutations = parse_vcf_files(str(vcf_dir))

    if mutations.empty:
        raise RuntimeError(f"No mutations parsed from {vcf_dir}")

    print(f"Parsed {len(mutations)} mutation rows for {mutations['sample'].nunique()} samples.")

    # Optionally save the long-format mutation table
    mutations_path = outdir / "mutations_long.tsv"
    mutations.to_csv(mutations_path, sep="\t", index=False)
    print(f"Wrote mutations table to: {mutations_path}")

    # ------------------------------------------------------------------
    # 2) Build mutation count matrix (samples x STR features)
    # ------------------------------------------------------------------
    print("[2/7] Building mutation count matrix...")
    matrix = build_mutation_matrix(
        mutations,
        ru="length",      # use motif length: LEN1, LEN2, ...
        ref_length=True,  # include reference repeat length
        change=True,      # include tumor–normal change (somatic events only)
    )

    if matrix.empty:
        raise RuntimeError("Mutation matrix is empty after tallying.")

    print(f"Matrix shape: {matrix.shape[0]} samples x {matrix.shape[1]} features")

    matrix_path = outdir / "matrix_raw.tsv"
    matrix.to_csv(matrix_path, sep="\t")
    print(f"Wrote raw matrix to: {matrix_path}")

    # ------------------------------------------------------------------
    # 3) Filter matrix
    # ------------------------------------------------------------------
    print("[3/7] Filtering matrix (manual thresholds)...")
    matrix_filt, summary = filter_mutation_matrix(
        matrix,
        feature_method="manual",
        min_feature_total=1,
        min_samples_with_feature=1,
        min_sample_total=1,
    )

    if matrix_filt.empty:
        raise RuntimeError("Filtered matrix is empty; filtering removed all rows/features.")

    print(
        f"Filtered matrix shape: {matrix_filt.shape[0]} samples x {matrix_filt.shape[1]} features\n"
        f"Original shape: {matrix.shape[0]} x {matrix.shape[1]}"
    )

    matrix_filt_path = outdir / "matrix_filtered.tsv"
    matrix_filt.to_csv(matrix_filt_path, sep="\t")
    print(f"Wrote filtered matrix to: {matrix_filt_path}")

    # ------------------------------------------------------------------
    # 4) Run NMF
    # ------------------------------------------------------------------
    print(f"[4/7] Running NMF with K={n_signatures} signatures...")
    nmf_res = run_nmf(
        matrix_filt,
        n_signatures=n_signatures,
        init="nndsvd",
        max_iter=200,
        random_state=0,
    )

    print("NMF finished.")
    print(f"Signatures shape: {nmf_res.signatures.shape}")
    print(f"Exposures shape:  {nmf_res.exposures.shape}")

    # ------------------------------------------------------------------
    # 5) Save & reload NMF result
    # ------------------------------------------------------------------
    nmf_dir = outdir / "nmf_results"
    print(f"[5/7] Saving NMF result to: {nmf_dir}")
    save_nmf_result(nmf_res, nmf_dir)

    nmf_loaded = load_nmf_result(nmf_dir)
    print("Reloaded NMF result from disk.")

    # ------------------------------------------------------------------
    # 6) Plot signatures, exposures, PCA
    # ------------------------------------------------------------------
    print("[6/7] Creating plots (signatures, exposures, PCA)...")

    fig_sig = plot_signatures(nmf_loaded, top_n=10)
    fig_sig_path = outdir / "signatures_top10.png"
    fig_sig.savefig(fig_sig_path, bbox_inches="tight")
    plt.close(fig_sig)

    fig_exp = plot_exposures(nmf_loaded)
    fig_exp_path = outdir / "exposures.png"
    fig_exp.savefig(fig_exp_path, bbox_inches="tight")
    plt.close(fig_exp)

    coords, var_ratio = compute_pca(nmf_loaded.exposures, n_components=2)
    print(f"PCA explained variance ratio: {var_ratio}")

    ax_pca = plot_pca_samples(coords, title="PCA of NMF exposures")
    fig_pca = ax_pca.figure
    fig_pca_path = outdir / "pca_exposures.png"
    fig_pca.savefig(fig_pca_path, bbox_inches="tight")
    plt.close(fig_pca)

    print(f"Saved plots to:\n  {fig_sig_path}\n  {fig_exp_path}\n  {fig_pca_path}")

    # ------------------------------------------------------------------
    # 7) Project a *new* VCF onto learned signatures
    # ------------------------------------------------------------------
    print("[7/7] Projecting a new VCF onto learned signatures...")

    vcf_files = sorted(vcf_dir.glob("*.vcf*"))
    if not vcf_files:
        raise RuntimeError(f"No VCF files found in {vcf_dir}")

    first_vcf = vcf_files[0]
    print(f"Using VCF for projection example: {first_vcf.name}")

    # Process a single VCF file to simulate a new cohort
    rows = process_vcf_to_rows(str(first_vcf))
    new_mut_df = pd.DataFrame(rows)
    if new_mut_df.empty:
        raise RuntimeError(f"No STR rows extracted from {first_vcf}")

    new_matrix = build_mutation_matrix(
        new_mut_df,
        ru="length",
        ref_length=True,
        change=True,
    )
    if new_matrix.empty:
        raise RuntimeError("New matrix for projection is empty.")

    exposures_new = project_onto_signatures(
        new_matrix=new_matrix,
        signatures=nmf_loaded.signatures,
        method="nnls",
    )

    exposures_new_path = outdir / "new_exposures_single_vcf.tsv"
    exposures_new.to_csv(exposures_new_path, sep="\t")
    print(f"Wrote projected exposures for new sample(s) to: {exposures_new_path}")

    print("Done.")


if __name__ == "__main__":
    # Adjust these paths to your environment before running
    example_vcf_dir = "/path/to/your/directory"  # directory with STR-annotated paired tumor–normal VCFs
    example_outdir = "example_output"

    run_full_pipeline(example_vcf_dir, example_outdir, n_signatures=2)
