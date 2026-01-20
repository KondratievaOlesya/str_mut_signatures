"""Console script for str_mut_signatures."""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

import pandas as pd

from . import __version__
from .extract_tally.extract_mutations import parse_vcf_files
from .extract_tally.filter import filter_mutation_matrix
from .extract_tally.tally import build_mutation_matrix
from .nmf.nmf import (
    load_nmf_result,
    project_onto_signatures,
    save_nmf_result,
)
from .nmf.nmf import (
    run_nmf as run_nmf_model,
)


class ValidationError(Exception):
    """Custom validation error for CLI argument checks."""

    pass


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for CLI."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def create_parser() -> argparse.ArgumentParser:
    """Create CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="STR mutation signature analysis from paired tumor–normal VCF files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # 1) Extract somatic STR mutation counts from paired tumor–normal VCFs
  str_mut_signatures extract \
    --vcf-dir data/vcfs \
    --out-matrix counts_len.tsv \
    --ru-length \
    --ref-length \
    --change

  # 2) Filter a count matrix (feature-level filtering)
  str_mut_signatures filter \\
      --matrix counts_len1.tsv \\
      --out-matrix counts_len1.filtered.tsv \\
      --feature-method elbow

  # 3) Run NMF decomposition on a filtered matrix and save signatures/exposures
  str_mut_signatures nmf \\
      --matrix counts_len1.filtered.tsv \\
      --outdir nmf_results \\
      --n-signatures 5

  # 4) Project a new cohort onto pre-computed signatures
  str_mut_signatures project \\
      --matrix new_counts.tsv \\
      --nmf-dir nmf_results \\
      --out-exposures new_exposures.tsv

  # Enable verbose logging
  str_mut_signatures extract \\
      --vcf-dir data/vcfs/ \\
      --out-matrix counts_len1.tsv \\
      --ru-length --ref-length --change \\
      --verbose
        """,
    )

    # Global options
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
        help="Subcommands",
    )

    # ------------------------------------------------------------------
    # extract subcommand
    # ------------------------------------------------------------------
    extract_parser = subparsers.add_parser(
        "extract",
        help="Extract somatic STR mutation counts from VCFs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Extract somatic STR mutation counts from paired tumor–normal VCF files.",
    )

    extract_parser.add_argument(
        "--vcf-dir",
        required=True,
        type=str,
        help="Directory with STR-annotated, paired tumor–normal VCF files.",
    )

    extract_parser.add_argument(
        "--out-matrix",
        required=True,
        type=str,
        help="Path to output TSV file with samples as rows and STR mutation features as columns.",
    )

    extract_parser.add_argument(
        "--ru-length",
        action="store_true",
        default=False,
        help="Include repeat-unit length as LEN{len(motif)} in feature labels.",
    )

    extract_parser.add_argument(
        "--ru",
        choices=["class", "ru"],
        default=None,
        help=(
            "How to include repeat-unit content in feature labels: "
            "'class' (base class AT/GC/MX) or "
            "'ru' (full repeat-unit sequence). "
            "If not specified, repeat-unit content is not included."
        ),
    )

    extract_parser.add_argument(
        "--ref-length",
        action="store_true",
        help="Include reference repeat length in feature labels.",
    )

    extract_parser.add_argument(
        "--change",
        action="store_true",
        help="Encode tumor–normal repeat-length change and restrict to somatic events.",
    )

    # ------------------------------------------------------------------
    # filter subcommand
    # ------------------------------------------------------------------
    filter_parser = subparsers.add_parser(
        "filter",
        help="Filter a STR mutation count matrix (features / samples).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Filter a STR mutation count matrix using different heuristics.",
    )

    filter_parser.add_argument(
        "--matrix",
        required=True,
        type=str,
        help="Input TSV count matrix (samples x STR mutation features).",
    )

    filter_parser.add_argument(
        "--out-matrix",
        required=True,
        type=str,
        help="Output TSV path for the filtered matrix.",
    )

    filter_parser.add_argument(
        "--feature-method",
        choices=["manual", "elbow", "percentile"],
        default="manual",
        help=(
            "Feature filtering method: "
            "'manual' (use explicit thresholds), "
            "'elbow' (elbow heuristic on feature totals), "
            "'percentile' (keep features above a percentile of totals). "
            "Default: manual"
        ),
    )

    filter_parser.add_argument(
        "--min-feature-total",
        type=int,
        default=10,
        help="Minimum total count across all samples for a feature (manual mode).",
    )

    filter_parser.add_argument(
        "--min-samples-with-feature",
        type=int,
        default=3,
        help="Minimum number of samples in which the feature must be non-zero.",
    )

    filter_parser.add_argument(
        "--min-sample-total",
        type=int,
        default=0,
        help="Minimum total count per sample (rows with less are dropped).",
    )

    filter_parser.add_argument(
        "--feature-percentile",
        type=float,
        default=0.9,
        help=(
            "Percentile (0–1) of feature totals used as threshold in "
            "feature-method=percentile. Default: 0.9"
        ),
    )

    # ------------------------------------------------------------------
    # nmf subcommand
    # ------------------------------------------------------------------
    nmf_parser = subparsers.add_parser(
        "nmf",
        help="Run NMF on a STR mutation count matrix",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Run NMF-based STR mutation signature decomposition on a count matrix.",
    )

    nmf_parser.add_argument(
        "--matrix",
        required=True,
        type=str,
        help="Input TSV count matrix (samples x STR mutation features).",
    )

    nmf_parser.add_argument(
        "--outdir",
        required=True,
        type=str,
        help="Output directory for NMF results (signatures.tsv, exposures.tsv, metadata.json).",
    )

    nmf_parser.add_argument(
        "--n-signatures",
        required=True,
        type=int,
        help="Number of signatures (rank) for NMF.",
    )

    nmf_parser.add_argument(
        "--max-iter",
        type=int,
        default=200,
        help="Maximum number of NMF iterations. Default: 200",
    )

    nmf_parser.add_argument(
        "--random-state",
        type=int,
        default=0,
        help="Random seed for NMF. Default: 0",
    )

    nmf_parser.add_argument(
        "--init",
        type=str,
        default="nndsvd",
        help="Initialization method for NMF (passed to sklearn.decomposition.NMF).",
    )

    nmf_parser.add_argument(
        "--alpha-W",
        type=float,
        default=0.0,
        help="L1/L2 regularization parameter for the W (exposure) matrix.",
    )

    nmf_parser.add_argument(
        "--alpha-H",
        type=float,
        default=0.0,
        help="L1/L2 regularization parameter for the H (signature) matrix.",
    )

    nmf_parser.add_argument(
        "--l1-ratio",
        type=float,
        default=0.0,
        help="Elastic-net mixing parameter (0 = L2, 1 = L1).",
    )

    # ------------------------------------------------------------------
    # project subcommand
    # ------------------------------------------------------------------
    project_parser = subparsers.add_parser(
        "project",
        help="Project new samples onto existing signatures (NMF result).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Given a new count matrix and an existing NMF result directory, "
            "compute exposures of new samples to the learned signatures."
        ),
    )

    project_parser.add_argument(
        "--matrix",
        required=True,
        type=str,
        help="Input TSV count matrix for NEW samples (samples x features).",
    )

    project_parser.add_argument(
        "--nmf-dir",
        required=True,
        type=str,
        help="Directory with a saved NMF result (signatures.tsv, metadata.json, ...).",
    )

    project_parser.add_argument(
        "--out-exposures",
        required=True,
        type=str,
        help="Output TSV for new sample exposures to the signatures.",
    )

    return parser


def validate_args(args: argparse.Namespace) -> None:
    """Validate CLI arguments."""
    if args.command == "extract":
        if not os.path.isdir(args.vcf_dir):
            raise ValidationError(f"--vcf-dir not found or not a directory: {args.vcf_dir}")

        out_parent = Path(args.out_matrix).resolve().parent
        if not out_parent.exists():
            raise ValidationError(f"Parent directory for --out-matrix not found: {out_parent}")

    elif args.command == "filter":
        if not os.path.isfile(args.matrix):
            raise ValidationError(f"--matrix file not found: {args.matrix}")

        out_parent = Path(args.out_matrix).resolve().parent
        if not out_parent.exists():
            raise ValidationError(f"Parent directory for --out-matrix not found: {out_parent}")

        if args.feature_method == "percentile" and not (0.0 <= args.feature_percentile <= 1.0):
            raise ValidationError("--feature-percentile must be between 0 and 1.")

    elif args.command == "nmf":
        if not os.path.isfile(args.matrix):
            raise ValidationError(f"--matrix file not found: {args.matrix}")

        if args.n_signatures <= 0:
            raise ValidationError("--n-signatures must be a positive integer")

        outdir = Path(args.outdir)
        out_parent = outdir.resolve().parent
        if not out_parent.exists():
            raise ValidationError(f"Parent directory for --outdir not found: {out_parent}")

    elif args.command == "project":
        if not os.path.isfile(args.matrix):
            raise ValidationError(f"--matrix file not found: {args.matrix}")

        nmf_dir = Path(args.nmf_dir)
        if not nmf_dir.is_dir():
            raise ValidationError(f"--nmf-dir not found or not a directory: {nmf_dir}")

        if not (nmf_dir / "signatures.tsv").is_file():
            raise ValidationError(
                f"--nmf-dir does not look like a saved NMF result "
                f"(missing signatures.tsv): {nmf_dir}"
            )

        out_parent = Path(args.out_exposures).resolve().parent
        if not out_parent.exists():
            raise ValidationError(f"Parent directory for --out-exposures not found: {out_parent}")


# ----------------------------------------------------------------------
# Command runners
# ----------------------------------------------------------------------


def run_extract_cli(args: argparse.Namespace, logger: logging.Logger) -> None:
    """Run the extract workflow."""
    logger.info("Loading VCFs from directory: %s", args.vcf_dir)
    mutations = parse_vcf_files(args.vcf_dir)

    if mutations.empty:
        raise RuntimeError(f"No mutations parsed from directory: {args.vcf_dir}")

    logger.info("Building mutation count matrix...")

    matrix = build_mutation_matrix(
        mutations,
        ru_length=args.ru_length,
        ru=args.ru,  # None | "class" | "ru"
        ref_length=args.ref_length,
        change=args.change,
    )

    if matrix.empty:
        raise RuntimeError("Resulting mutation matrix is empty. Check input VCFs and options.")

    logger.info("Writing count matrix to: %s", args.out_matrix)
    matrix.to_csv(args.out_matrix, sep="\t")
    logger.info("Matrix shape: %s samples x %s features", matrix.shape[0], matrix.shape[1])


def run_filter_cli(args: argparse.Namespace, logger: logging.Logger) -> None:
    """Run the filter workflow."""
    logger.info("Loading count matrix from: %s", args.matrix)
    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)

    if matrix.empty:
        raise RuntimeError("Input matrix is empty, cannot filter.")

    logger.info(
        "Filtering matrix with method=%s (min_feature_total=%d, "
        "min_samples_with_feature=%d, min_sample_total=%d, feature_percentile=%.3f)",
        args.feature_method,
        args.min_feature_total,
        args.min_samples_with_feature,
        args.min_sample_total,
        args.feature_percentile,
    )

    filtered, summary = filter_mutation_matrix(
        matrix,
        feature_method=args.feature_method,
        min_feature_total=args.min_feature_total,
        min_samples_with_feature=args.min_samples_with_feature,
        min_sample_total=args.min_sample_total,
        feature_percentile=args.feature_percentile,
    )

    if filtered.empty:
        raise RuntimeError("Filtered matrix is empty; filtering removed all rows/features.")

    # Use shapes; don't assume particular summary attribute names
    logger.info(
        "Filtered matrix shape: %d samples x %d features (original: %d samples x %d features).",
        filtered.shape[0],
        filtered.shape[1],
        matrix.shape[0],
        matrix.shape[1],
    )

    logger.info("Writing filtered matrix to: %s", args.out_matrix)
    filtered.to_csv(args.out_matrix, sep="\t")


def run_nmf_cli(args: argparse.Namespace, logger: logging.Logger) -> None:
    """Run the NMF workflow."""
    logger.info("Loading count matrix from: %s", args.matrix)
    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)

    if matrix.empty:
        raise RuntimeError("Input matrix is empty, cannot run NMF.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Running NMF with K=%d (init=%s, max_iter=%d, random_state=%d, alpha_W=%.3f, alpha_H=%.3f, l1_ratio=%.3f)...",
        args.n_signatures,
        args.init,
        args.max_iter,
        args.random_state,
        args.alpha_W,
        args.alpha_H,
        args.l1_ratio,
    )

    nmf_res = run_nmf_model(
        matrix,
        n_signatures=args.n_signatures,
        init=args.init,
        max_iter=args.max_iter,
        random_state=args.random_state,
        alpha_W=args.alpha_W,
        alpha_H=args.alpha_H,
        l1_ratio=args.l1_ratio,
    )

    # Save signatures/exposures/metadata in a versioned, stable format
    save_nmf_result(nmf_res, outdir)

    logger.info(
        "NMF completed. Signatures: %s features x %s signatures; exposures: %s samples x %s signatures.",
        nmf_res.signatures.shape[0],
        nmf_res.signatures.shape[1],
        nmf_res.exposures.shape[0],
        nmf_res.exposures.shape[1],
    )

    logger.info("Saved NMF result to: %s", outdir)


def run_project_cli(args: argparse.Namespace, logger: logging.Logger) -> None:
    """Run the projection workflow (new samples -> existing signatures)."""
    logger.info("Loading new count matrix from: %s", args.matrix)
    new_matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)

    if new_matrix.empty:
        raise RuntimeError("Input matrix for projection is empty.")

    nmf_dir = Path(args.nmf_dir)
    logger.info("Loading NMF result from: %s", nmf_dir)
    nmf_res = load_nmf_result(nmf_dir)

    logger.info(
        "Projecting %d new samples onto %d signatures...",
        new_matrix.shape[0],
        nmf_res.signatures.shape[1],
    )

    exposures_new = project_onto_signatures(
        new_matrix=new_matrix,
        signatures=nmf_res.signatures,
        method="nnls",
    )

    logger.info("Writing new exposures to: %s", args.out_exposures)
    exposures_new.to_csv(args.out_exposures, sep="\t")


# ----------------------------------------------------------------------
# main entry
# ----------------------------------------------------------------------


def main() -> int:
    """CLI entry point with argument parsing and validation."""
    parser = create_parser()
    args = parser.parse_args()

    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    try:
        # Validate arguments
        validate_args(args)

        if args.command == "extract":
            run_extract_cli(args, logger)
        elif args.command == "filter":
            run_filter_cli(args, logger)
        elif args.command == "nmf":
            run_nmf_cli(args, logger)
        elif args.command == "project":
            run_project_cli(args, logger)
        else:
            raise ValidationError(f"Unknown command: {args.command}")

        logger.info("Done.")
        return 0

    except ValidationError as e:
        logger.error("Validation error: %s", e)
        return 1

    except FileNotFoundError as e:
        logger.error("File not found: %s", e)
        return 1

    except KeyboardInterrupt:
        logger.error("Interrupted by user.")
        return 1

    except Exception as e:
        # Show traceback only in verbose mode
        logger.error("Unexpected error: %s", e, exc_info=args.verbose)
        return 1


if __name__ == "__main__":
    sys.exit(main())
