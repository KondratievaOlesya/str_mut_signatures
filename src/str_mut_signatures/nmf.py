import os
import pandas as pd
import numpy as np


def run_nmf_decomposition(counts_matrix_file, output_dir, n_components=3, random_state=42):
    """
    Run NMF decomposition on counts matrix.

    Args:
        counts_matrix_file (str): Path to counts matrix TSV file
        output_dir (str): Output directory for results
        n_components (int): Number of components for NMF
        random_state (int): Random state for reproducibility
    """

    # Check if sklearn is available
    try:
        from sklearn.decomposition import NMF
    except ImportError as e:
        raise ImportError(
            "scikit-learn is required for NMF. Install with: pip install scikit-learn"
        ) from e

    # Check if plotting libraries are available
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        print("Warning: matplotlib or seaborn not available. Plots will be skipped.")
        PLOTTING_AVAILABLE = False
    else:
        PLOTTING_AVAILABLE = True

    os.makedirs(output_dir, exist_ok=True)

    # Read counts matrix
    counts_df = pd.read_csv(counts_matrix_file, sep='\t', index_col=0)

    if counts_df.empty:
        raise ValueError("Counts matrix is empty")

    # Remove any columns with all zeros
    counts_df = counts_df.loc[:, (counts_df != 0).any(axis=0)]

    if counts_df.empty:
        raise ValueError("No non-zero columns in counts matrix")

    print(f"Running NMF on matrix shape: {counts_df.shape}")
    print(f"Number of components: {n_components}")

    # Convert to numpy array
    X = counts_df.values

    # Run NMF
    nmf = NMF(n_components=n_components, random_state=random_state, max_iter=1000)
    W = nmf.fit_transform(X)  # Sample signatures
    H = nmf.components_  # Signature profiles

    # Create results DataFrames
    signatures_df = pd.DataFrame(
        H,
        columns=counts_df.columns,
        index=[f'Signature_{i + 1}' for i in range(n_components)]
    )

    exposures_df = pd.DataFrame(
        W,
        index=counts_df.index,
        columns=[f'Signature_{i + 1}' for i in range(n_components)]
    )

    # Save results
    signatures_df.to_csv(os.path.join(output_dir, 'signatures.tsv'), sep='\t')
    exposures_df.to_csv(os.path.join(output_dir, 'exposures.tsv'), sep='\t')

    # Save reconstruction error
    with open(os.path.join(output_dir, 'nmf_metrics.txt'), 'w') as f:
        f.write(f"Reconstruction error: {nmf.reconstruction_err_:.6f}\n")
        f.write(f"Number of components: {n_components}\n")
        f.write(f"Input matrix shape: {counts_df.shape}\n")

    # Create plots if available
    if PLOTTING_AVAILABLE:
        try:
            plt.figure(figsize=(12, 8))

            # Plot signatures
            plt.subplot(2, 1, 1)
            signatures_df.T.plot(kind='bar', ax=plt.gca())
            plt.title('Mutation Signatures')
            plt.ylabel('Contribution')
            plt.legend(title='Signatures')
            plt.xticks(rotation=45)

            # Plot exposures
            plt.subplot(2, 1, 2)
            exposures_df.plot(kind='bar', stacked=True, ax=plt.gca())
            plt.title('Sample Exposures to Signatures')
            plt.ylabel('Exposure')
            plt.xlabel('Samples')
            plt.legend(title='Signatures')
            plt.xticks(rotation=45)

            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'nmf_results.png'),
                        dpi=300, bbox_inches='tight')
            plt.close()

            print("Plots saved as nmf_results.png")

        except Exception as e:
            print(f"Warning: Could not create plots: {e}")

    print(f"NMF analysis completed. Results saved to: {output_dir}")
    print(f"Reconstruction error: {nmf.reconstruction_err_:.6f}")


def run_sigprofiler_extractor(counts_matrix_file, output_dir):
    """
    Run SigProfilerExtractor if available.
    """
    try:
        from SigProfilerExtractor import sigpro as sig
    except ImportError:
        print("SigProfilerExtractor not available. Using sklearn NMF instead.")
        return run_nmf_decomposition(counts_matrix_file, output_dir)

    # Read and prepare data
    counts_df = pd.read_csv(counts_matrix_file, sep='\t', index_col=0)
    counts_t = counts_df.T

    # Save temporary input
    temp_input = os.path.join(output_dir, "sigprofiler_input.txt")
    counts_t.to_csv(temp_input, sep='\t')

    print("Running SigProfilerExtractor...")

    # Run SigProfiler with minimal parameters for testing
    sig.sigProfilerExtractor(
        "matrix",
        output_dir,
        temp_input,
        reference_genome="GRCh37",
        context_type="custom",
        minimum_signatures=1,
        maximum_signatures=5,
        nmf_replicates=5,
        cpu=-1
    )

    print(f"SigProfilerExtractor completed. Results in: {output_dir}")