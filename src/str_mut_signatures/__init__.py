"""
Top-level package for str_mut_signatures.

Analysis of Short Tandem Repeat (STR) mutation signatures from VCF files.
Provides both library and CLI interfaces.
"""

from __future__ import annotations

from .extract_tally.extract_mutations import (
    parse_vcf_files,
    process_vcf_to_rows,
    save_counts_matrix,
)
from .extract_tally.filter import (
    filter_mutation_matrix,
)
from .extract_tally.tally import (
    build_mutation_matrix,
)
from .nmf.nmf import (
    NMFResult,
    load_nmf_result,
    project_onto_signatures,
    run_nmf,
    save_nmf_result,
)
from .nmf.plot import (
    compute_pca,
    plot_exposures,
    plot_pca_samples,
    plot_signatures,
)

__author__ = "Olesia Kondrateva"
__email__ = "xkdnoa@gmail.com"
__version__ = "2.0.0"

__all__ = [
    # VCF → mutations
    "parse_vcf_files",
    "save_counts_matrix",
    "process_vcf_to_rows",
    # mutations → matrix
    "build_mutation_matrix",
    # matrix filtering
    "filter_mutation_matrix",
    # NMF core
    "NMFResult",
    "run_nmf",
    # NMF persistence
    "save_nmf_result",
    "load_nmf_result",
    # applying signatures to new data
    "project_onto_signatures",
    # Plotting
    "compute_pca",
    "plot_exposures",
    "plot_pca_samples",
    "plot_signatures"
]
