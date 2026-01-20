from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from sklearn.cluster import KMeans
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score

try:
    # Python 3.11+
    from datetime import UTC
except ImportError:  # Python < 3.11
    from datetime import timezone

    UTC = timezone.utc
FORMAT_VERSION = 1


@dataclass
class NMFResult:
    """
    Container for NMF-based STR mutation signature decomposition.

    Attributes
    ----------
    signatures : pandas.DataFrame
        :no-index:
        Matrix of signature profiles.

        - index   : features (same as input ``matrix.columns``)
        - columns : signatures (``Signature_1``, ``Signature_2``, ...)

    exposures : pandas.DataFrame
        :no-index:
        Matrix of sample exposures to each signature.

        - index   : samples (same as input ``matrix.index``)
        - columns : signatures (``Signature_1``, ``Signature_2``, ...)

    groups : pandas.DataFrame
        :no-index:
        Sample-level grouping or annotation table aligned to exposures.
        Typically indexed by sample (same as input ``matrix.index``).

    model_params : dict[str, Any]
        :no-index:
        Hyperparameters and metadata used to fit the model
        (e.g. ``n_signatures``, ``init``, ``max_iter``, ``random_state``).
    """

    signatures: pd.DataFrame
    exposures: pd.DataFrame
    groups: pd.DataFrame
    model_params: dict[str, Any]


def make_groups_df(
    sample_index: pd.Index,
    cluster_labels: np.ndarray | None = None,
) -> pd.DataFrame:
    """
    Create a groups DataFrame indexed by samples.

    The returned DataFrame always contains at least one column:

    - ``group``: string labels for each sample (default: ``"1"`` for all samples)

    If ``cluster_labels`` is provided, the ``group`` column is overwritten with
    the corresponding cluster identifiers (converted to strings).

    Parameters
    ----------
    sample_index : pandas.Index
        Index of sample identifiers to use for the resulting DataFrame.
    cluster_labels : numpy.ndarray or None, optional
        Array of cluster labels aligned with ``sample_index``.
        If provided, these labels are used to populate the ``group`` column.

    Returns
    -------
    pandas.DataFrame
        Groups DataFrame indexed by sample, with at least a ``"group"`` column.
    """
    groups = pd.DataFrame(index=sample_index)
    groups["group"] = "1"

    if cluster_labels is not None:
        if len(cluster_labels) != len(sample_index):
            raise ValueError(
                "cluster_labels length must match number of samples "
                f"({len(cluster_labels)} != {len(sample_index)})."
            )
        groups["group"] = pd.Series(cluster_labels, index=sample_index).astype(str)

    return groups


def add_labels_to_groups(
    groups_df: pd.DataFrame,
    labels: pd.Series | pd.DataFrame,
    *,
    label_col: str = "label",
    allow_missing: bool = True,
) -> pd.DataFrame:
    """
    Add external labels to a groups DataFrame by sample index.

    Labels are merged into ``groups_df`` using the sample identifier as index.
    All samples present in ``groups_df`` are preserved in the output.

    The ``labels`` input may be provided in one of the following forms:

    - :class:`pandas.Series` with sample IDs as index.
    - :class:`pandas.DataFrame` with sample IDs as index and either:
        * a single column, or
        * a column named according to ``label_col``.

    Parameters
    ----------
    groups_df : pandas.DataFrame
        Groups DataFrame indexed by sample identifier.
    labels : pandas.Series or pandas.DataFrame
        External labels aligned by sample index.
    label_col : str, optional
        Name of the column to use or create for the merged labels.
        Default is ``"label"``.
    allow_missing : bool, optional
        Whether to allow missing labels for some samples.
        If ``False``, a :class:`ValueError` is raised when missing labels are found.

    Returns
    -------
    pandas.DataFrame
        Copy of ``groups_df`` with an additional ``label_col`` column
        containing the merged labels.

    Raises
    ------
    ValueError
        If ``allow_missing`` is ``False`` and some samples in ``groups_df``
        have no corresponding label in ``labels``.
    """
    out = groups_df.copy()

    if isinstance(labels, pd.Series):
        lab = labels.rename(label_col)
    elif isinstance(labels, pd.DataFrame):
        if label_col in labels.columns:
            lab = labels[label_col]
        elif labels.shape[1] == 1:
            lab = labels.iloc[:, 0].rename(label_col)
        else:
            raise ValueError(
                f"labels DataFrame must have column '{label_col}' or exactly one column."
            )
    else:
        raise TypeError("labels must be a pandas Series or DataFrame.")

    # Align by index (sample id)
    out[label_col] = lab.reindex(out.index)

    if not allow_missing and out[label_col].isna().any():
        missing = out.index[out[label_col].isna()]
        raise ValueError(
            f"Missing labels for {len(missing)} samples (first 10): {missing[:10].tolist()}"
        )

    return out


def cluster_samples(
    exposures: pd.DataFrame,
    max_clusters: int = 6,
    random_state: int | None = 0,
) -> np.ndarray | None:
    """
    Cluster samples using KMeans and select the best ``k`` via silhouette score.

    Clustering is performed on the rows of ``exposures`` (one row per sample).
    The number of clusters is selected by evaluating silhouette scores for
    candidate values of ``k`` up to ``max_clusters``.

    If clustering cannot be performed (e.g., too few samples or invalid input),
    the function returns ``None``.

    Parameters
    ----------
    exposures : pandas.DataFrame
        Sample-by-signature exposure matrix. Rows correspond to samples and
        columns correspond to signatures.
    max_clusters : int, optional
        Maximum number of clusters to consider. Default is 6.
    random_state : int or None, optional
        Random seed passed to KMeans for reproducibility. If ``None``, the
        estimator is not seeded.

    Returns
    -------
    numpy.ndarray or None
        Array of cluster labels (``0`` .. ``k-1``) in the same order as
        ``exposures.index``, or ``None`` if clustering is not possible.
    """
    n_samples = exposures.shape[0]
    if n_samples < 3:
        return None

    X = exposures.to_numpy(dtype=float)
    best_k = None
    best_score = -np.inf
    best_labels = None

    for k in range(2, min(max_clusters, n_samples - 1) + 1):
        km = KMeans(n_clusters=k, n_init="auto", random_state=random_state)
        labels = km.fit_predict(X)

        if len(np.unique(labels)) < 2:
            continue

        try:
            score = silhouette_score(X, labels)
        except ValueError:
            continue

        if score > best_score:
            best_score = score
            best_k = k
            best_labels = labels

    return best_labels if best_k is not None else None


def validate_input_matrix(matrix: pd.DataFrame) -> np.ndarray:
    """
    Validate that the input is a non-empty, numeric, non-negative matrix.

    This function checks that ``matrix`` is:

    - A :class:`pandas.DataFrame`
    - Non-empty
    - Fully numeric
    - Contains no negative values

    On success, the underlying numeric values are returned as a NumPy array.

    Parameters
    ----------
    matrix : pandas.DataFrame
        Input matrix with samples as rows and features as columns.

    Returns
    -------
    numpy.ndarray
        The underlying numeric array with shape
        ``(n_samples, n_features)``.

    Raises
    ------
    ValueError
        If the input is empty, contains non-numeric values, or contains
        negative entries.
    """
    if not isinstance(matrix, pd.DataFrame):
        raise TypeError("matrix must be a pandas.DataFrame")

    if matrix.empty:
        raise ValueError("matrix is empty; NMF requires a non-empty matrix.")

    # Ensure numeric dtype
    if not np.issubdtype(matrix.dtypes.values[0], np.number) or not all(
        np.issubdtype(dtype, np.number) for dtype in matrix.dtypes.values
    ):
        raise TypeError("matrix must contain only numeric values.")

    values = matrix.to_numpy(dtype=float)

    if (values < 0).any():
        raise ValueError("NMF requires a non-negative matrix; found negative entries.")

    return values


def run_nmf(
    matrix: pd.DataFrame,
    n_signatures: int,
    init: str = "nndsvd",
    max_iter: int = 200,
    random_state: int | None = 0,
    alpha_W: float = 0.0,
    alpha_H: float = 0.0,
    l1_ratio: float = 0.0,
    max_clusters: int = 1,  # <=1 means "no clustering"
) -> NMFResult:
    """
    Run NMF decomposition on an STR mutation count matrix.

    This function factorizes a non-negative mutation count matrix into:

    - signature profiles (feature-by-signature)
    - sample exposures (sample-by-signature)

    Optionally, samples can be clustered based on their exposure profiles.

    Parameters
    ----------
    matrix : pandas.DataFrame
        Non-negative count matrix.

        - rows    : samples
        - columns : mutation feature categories
    n_signatures : int
        Number of signatures (components) to extract.
    init : str, optional
        Initialization method for NMF (passed to the underlying estimator).
        Default is ``"nndsvd"``.
    max_iter : int, optional
        Maximum number of iterations. Default is 200.
    random_state : int or None, optional
        Random seed for reproducibility. If ``None``, the estimator is not seeded.
        Default is 0.
    alpha_W : float, optional
        Regularization parameter for the W matrix (exposures), if supported by
        the chosen NMF implementation. Default is 0.0.
    alpha_H : float, optional
        Regularization parameter for the H matrix (signatures), if supported by
        the chosen NMF implementation. Default is 0.0.
    l1_ratio : float, optional
        The Elastic-Net mixing parameter, with ``0 <= l1_ratio <= 1``.
        Default is 0.0.
    max_clusters : int, optional
        Maximum number of clusters to consider for optional exposure-based
        clustering. Values ``<= 1`` disable clustering. Default is 1.

    Returns
    -------
    NMFResult
        Container with signature profiles, exposures, optional grouping
        information, and model parameters.

    Raises
    ------
    ValueError
        If ``matrix`` is empty, contains non-numeric values, or contains
        negative entries.
    ValueError
        If ``n_signatures`` is not a positive integer.
    """
    values = validate_input_matrix(matrix)

    n_samples, n_features = values.shape
    if n_signatures <= 0:
        raise ValueError("n_signatures must be a positive integer.")
    if n_signatures > min(n_samples, n_features):
        raise ValueError(
            f"n_signatures ({n_signatures}) cannot exceed "
            f"min(n_samples, n_features) = {min(n_samples, n_features)}."
        )

    model = NMF(
        n_components=n_signatures,
        init=init,
        max_iter=max_iter,
        random_state=random_state,
        alpha_W=alpha_W,
        alpha_H=alpha_H,
        l1_ratio=l1_ratio,
    )

    # W: samples x K (exposures), H: K x features (signatures)
    W = model.fit_transform(values)
    H = model.components_

    signature_labels = [f"Signature_{k + 1}" for k in range(n_signatures)]

    signatures_df = pd.DataFrame(
        H.T,
        index=matrix.columns,
        columns=signature_labels,
    )

    exposures_df = pd.DataFrame(
        W,
        index=matrix.index,
        columns=signature_labels,
    )

    cluster_labels = None
    if max_clusters is not None and max_clusters > 1:
        cluster_labels = cluster_samples(
            exposures_df,
            max_clusters=max_clusters,
            random_state=random_state,
        )
    groups_df = make_groups_df(exposures_df.index, cluster_labels=cluster_labels)

    model_params: dict[str, Any] = {
        "n_signatures": n_signatures,
        "init": init,
        "max_iter": max_iter,
        "random_state": random_state,
        "alpha_W": alpha_W,
        "alpha_H": alpha_H,
        "l1_ratio": l1_ratio,
        "reconstruction_err_": float(getattr(model, "reconstruction_err_", np.nan)),
        "n_iter_": int(getattr(model, "n_iter_", -1)),
        "format_version": FORMAT_VERSION,
        "max_clusters": int(max_clusters) if max_clusters is not None else None,
        "n_groups": int(groups_df["group"].nunique()),
        "created_at": datetime.now(UTC).isoformat(),
    }

    return NMFResult(
        signatures=signatures_df,
        exposures=exposures_df,
        groups=groups_df,
        model_params=model_params,
    )


# ---------------------------------------------------------------------------
# Saving / loading
# ---------------------------------------------------------------------------


def save_nmf_result(result: NMFResult, outdir: str | Path) -> None:
    """
    Save an :class:`NMFResult` to a directory on disk.

    This function writes the main components of the NMF decomposition to
    tabular and JSON files in the specified output directory.

    The following files are created:

    - ``signatures.tsv``:
        Signature profiles (features × K).
    - ``exposures.tsv``:
        Sample exposures (samples × K).
    - ``metadata.json``:
        JSON file containing ``model_params`` together with basic shape
        information and a ``format_version`` field.

    Parameters
    ----------
    result : NMFResult
        Result object containing signatures, exposures, groups, and model
        parameters to be saved.
    outdir : str or pathlib.Path
        Output directory where result files will be written.

    Returns
    -------
    None
    """
    outpath = Path(outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    sig_path = outpath / "signatures.tsv"
    exp_path = outpath / "exposures.tsv"
    meta_path = outpath / "metadata.json"
    groups_path = outpath / "groups.tsv"

    result.signatures.to_csv(sig_path, sep="\t")
    result.exposures.to_csv(exp_path, sep="\t")
    result.groups.to_csv(groups_path, sep="\t")

    # Add a bit of structural metadata
    metadata = dict(result.model_params)
    metadata.setdefault("format_version", FORMAT_VERSION)
    metadata["signatures_shape"] = list(result.signatures.shape)
    metadata["exposures_shape"] = list(result.exposures.shape)
    metadata["signature_columns"] = list(result.signatures.columns)
    metadata["exposure_columns"] = list(result.exposures.columns)
    metadata["groups_columns"] = list(result.groups.columns)

    meta_path.write_text(json.dumps(metadata, indent=2))


def load_nmf_result(outdir: str | Path) -> NMFResult:
    """
    Load an :class:`NMFResult` previously saved with :func:`save_nmf_result`.

    This function reads the files created by ``save_nmf_result`` from the given
    directory and reconstructs the corresponding :class:`NMFResult` object.

    The following files are expected in ``outdir``:

    - ``signatures.tsv``:
        Signature profiles (features × K).
    - ``exposures.tsv``:
        Sample exposures (samples × K).
    - ``metadata.json``:
        JSON file containing model parameters and basic shape information.

    Parameters
    ----------
    outdir : str or pathlib.Path
        Directory containing the saved NMF result files.

    Returns
    -------
    NMFResult
        Reconstructed NMF result with signatures, exposures, groups (if present),
        and model parameters.
    """
    outpath = Path(outdir)

    sig_path = outpath / "signatures.tsv"
    exp_path = outpath / "exposures.tsv"
    meta_path = outpath / "metadata.json"
    groups_path = outpath / "groups.tsv"

    if not sig_path.is_file():
        raise FileNotFoundError(f"Missing signatures file: {sig_path}")
    if not exp_path.is_file():
        raise FileNotFoundError(f"Missing exposures file: {exp_path}")
    if not meta_path.is_file():
        raise FileNotFoundError(f"Missing metadata file: {meta_path}")

    signatures = pd.read_csv(sig_path, sep="\t", index_col=0)
    exposures = pd.read_csv(exp_path, sep="\t", index_col=0)

    metadata = json.loads(meta_path.read_text())

    # Basic sanity checks (non-fatal: warn instead of raising)
    if "format_version" in metadata and metadata["format_version"] != FORMAT_VERSION:
        warnings.warn(
            f"NMF metadata format_version={metadata['format_version']} "
            f"differs from current FORMAT_VERSION={FORMAT_VERSION}. "
            "Results should still load, but fields may have changed.",
            RuntimeWarning,
            stacklevel=2,
        )

    # Groups: load if present; otherwise default to "1"
    if groups_path.is_file():
        groups = pd.read_csv(groups_path, sep="\t", index_col=0)
        # enforce alignment
        groups = groups.reindex(exposures.index)
        if "group" not in groups.columns:
            groups["group"] = "1"
        groups["group"] = groups["group"].fillna("1").astype(str)
    else:
        groups = make_groups_df(exposures.index)

    if "signatures_shape" in metadata:
        expected = tuple(metadata["signatures_shape"])
        if signatures.shape != expected:
            warnings.warn(
                f"Loaded signatures shape {signatures.shape} does not match metadata {expected}.",
                RuntimeWarning,
                stacklevel=2,
            )

    if "exposures_shape" in metadata:
        expected = tuple(metadata["exposures_shape"])
        if exposures.shape != expected:
            warnings.warn(
                f"Loaded exposures shape {exposures.shape} does not match metadata {expected}.",
                RuntimeWarning,
                stacklevel=2,
            )

    return NMFResult(
        signatures=signatures,
        exposures=exposures,
        groups=groups,
        model_params=metadata,
    )


# ---------------------------------------------------------------------------
# Projection of new data onto existing signatures
# ---------------------------------------------------------------------------
def validate_projection_inputs(
    new_matrix: pd.DataFrame,
    signatures: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    if not isinstance(new_matrix, pd.DataFrame):
        raise TypeError("new_matrix must be a pandas.DataFrame.")
    if not isinstance(signatures, pd.DataFrame):
        raise TypeError("signatures must be a pandas.DataFrame.")

    if new_matrix.empty:
        raise ValueError("new_matrix is empty; nothing to project.")
    if signatures.empty:
        raise ValueError("signatures is empty; nothing to project onto.")

    if not all(np.issubdtype(dtype, np.number) for dtype in new_matrix.dtypes.values):
        raise TypeError("new_matrix must contain only numeric values.")
    if not all(np.issubdtype(dtype, np.number) for dtype in signatures.dtypes.values):
        raise TypeError("signatures must contain only numeric values.")

    if (new_matrix.to_numpy(dtype=float) < 0).any():
        raise ValueError("new_matrix must be non-negative for NMF-based projection.")
    if (signatures.to_numpy(dtype=float) < 0).any():
        raise ValueError("signatures must be non-negative.")

    # Align features: new_matrix columns vs signatures index
    common_features = list(new_matrix.columns.intersection(signatures.index))
    if not common_features:
        raise ValueError("No overlapping features between new_matrix.columns and signatures.index.")

    # Warn if many are missing
    frac_missing = 1.0 - len(common_features) / len(signatures.index)
    if frac_missing > 0.5:
        warnings.warn(
            f"Only {len(common_features)} of {len(signatures.index)} signature features "
            f"are present in new_matrix (~{100 * (1 - frac_missing):.1f}% overlap).",
            RuntimeWarning,
            stacklevel=2,
        )

    # Subset and order both sides by common_features
    new_sub = new_matrix.loc[:, common_features]
    sig_sub = signatures.loc[common_features, :]

    return new_sub, sig_sub, common_features


def project_onto_signatures(
    new_matrix: pd.DataFrame,
    signatures: pd.DataFrame,
    method: str = "nnls",
) -> pd.DataFrame:
    """
    Project new samples onto existing signatures to obtain exposures.

    This function computes exposure weights for each new sample given a fixed
    set of signature profiles.

    Parameters
    ----------
    new_matrix : pandas.DataFrame
        Matrix of new samples.

        - rows    : samples
        - columns : features (must overlap ``signatures.index``)
    signatures : pandas.DataFrame
        Signature profile matrix.

        - index   : features (same feature space as ``new_matrix.columns``)
        - columns : signatures (e.g., ``"Signature_1"``, ``"Signature_2"``, ...)
    method : {"nnls"}, optional
        Projection method. Currently only non-negative least squares
        (``"nnls"``) is implemented.

    Returns
    -------
    pandas.DataFrame
        Exposure matrix for the new samples.

        - rows    : samples (same as ``new_matrix.index``)
        - columns : signatures (same as ``signatures.columns``)

    Notes
    -----
    For ``method="nnls"``, for each sample vector ``x`` (1 × F) the exposures
    ``e`` are obtained by solving::

        minimize || x - A e ||_2   subject to e >= 0

    where ``A`` is the feature-by-signature matrix (F × K).

    Raises
    ------
    ValueError
        If ``method`` is not supported or if there is no overlap between
        ``new_matrix.columns`` and ``signatures.index``.
    """
    if method != "nnls":
        raise ValueError('Only method="nnls" is currently supported.')

    if nnls is None:  # pragma: no cover
        raise ImportError(
            "scipy is required for method='nnls'. "
            "Install scipy or choose another projection method."
        )

    new_sub, sig_sub, common_features = validate_projection_inputs(new_matrix, signatures)

    # A: F x K
    A = sig_sub.to_numpy(dtype=float)
    K = A.shape[1]

    exposures = np.zeros((new_sub.shape[0], K), dtype=float)

    # Solve NNLS per sample
    for i, (_, row) in enumerate(new_sub.iterrows()):
        b = row.to_numpy(dtype=float)
        coeffs, _ = nnls(A, b)
        exposures[i, :] = coeffs

    exposures_df = pd.DataFrame(
        exposures,
        index=new_sub.index,
        columns=sig_sub.columns,
    )

    return exposures_df
