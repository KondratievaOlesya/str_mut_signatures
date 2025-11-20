from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from str_mut_signatures.nmf.nmf import (
    NMFResult,
    load_nmf_result,
    project_onto_signatures,
    run_nmf,
    save_nmf_result,
    validate_input_matrix,
)


class TestValidateInputMatrix:
    def test_non_dataframe_raises(self):
        with pytest.raises(TypeError):
            _ = validate_input_matrix("not a df")  # type: ignore[arg-type]

    def test_empty_dataframe_raises(self):
        df = pd.DataFrame()
        with pytest.raises(ValueError) as excinfo:
            _ = validate_input_matrix(df)
        assert "matrix is empty" in str(excinfo.value)

    def test_non_numeric_dataframe_raises(self):
        df = pd.DataFrame({"A": ["x", "y"], "B": ["1", "2"]})
        with pytest.raises(TypeError) as excinfo:
            _ = validate_input_matrix(df)
        assert "numeric" in str(excinfo.value)

    def test_negative_values_raise(self):
        df = pd.DataFrame({"A": [1, -1], "B": [0, 2]})
        with pytest.raises(ValueError) as excinfo:
            _ = validate_input_matrix(df)
        assert "non-negative matrix" in str(excinfo.value)

    def test_valid_matrix_returns_numpy_array(self):
        df = pd.DataFrame({"A": [1, 2], "B": [0, 3]})
        arr = validate_input_matrix(df)

        assert isinstance(arr, np.ndarray)
        assert arr.shape == (2, 2)
        assert arr.dtype == float
        assert np.all(arr >= 0)


class TestRunNMF:
    def toy_matrix(self) -> pd.DataFrame:
        # simple small non-negative matrix
        # samples x features
        return pd.DataFrame(
            {
                "f1": [1.0, 2.0, 3.0],
                "f2": [0.0, 1.0, 0.0],
                "f3": [4.0, 0.0, 1.0],
            },
            index=["s1", "s2", "s3"],
        )

    def test_basic_run_shapes_and_labels(self):
        matrix = self.toy_matrix()
        k = 2

        res = run_nmf(matrix, n_signatures=k, random_state=0)

        assert isinstance(res, NMFResult)

        # signatures: features x K
        assert isinstance(res.signatures, pd.DataFrame)
        assert res.signatures.shape == (matrix.shape[1], k)
        assert list(res.signatures.index) == list(matrix.columns)
        assert list(res.signatures.columns) == [f"Signature_{i+1}" for i in range(k)]

        # exposures: samples x K
        assert isinstance(res.exposures, pd.DataFrame)
        assert res.exposures.shape == (matrix.shape[0], k)
        assert list(res.exposures.index) == list(matrix.index)
        assert list(res.exposures.columns) == [f"Signature_{i+1}" for i in range(k)]

        # all non-negative (up to tiny numerical noise)
        assert np.all(res.signatures.to_numpy() >= -1e-10)
        assert np.all(res.exposures.to_numpy() >= -1e-10)

    def test_invalid_n_signatures_zero_or_negative(self):
        matrix = self.toy_matrix()

        with pytest.raises(ValueError):
            run_nmf(matrix, n_signatures=0)

        with pytest.raises(ValueError):
            run_nmf(matrix, n_signatures=-1)

    def test_invalid_n_signatures_too_large(self):
        matrix = self.toy_matrix()
        n_samples, n_features = matrix.shape
        too_large = min(n_samples, n_features) + 1

        with pytest.raises(ValueError) as excinfo:
            run_nmf(matrix, n_signatures=too_large)
        msg = str(excinfo.value)
        assert "cannot exceed min" in msg

    def test_reproducibility_with_random_state(self):
        matrix = self.toy_matrix()

        res1 = run_nmf(matrix, n_signatures=2, random_state=42)
        res2 = run_nmf(matrix, n_signatures=2, random_state=42)

        # Same results for same seed
        assert np.allclose(
            res1.signatures.to_numpy(), res2.signatures.to_numpy(), rtol=1e-6, atol=1e-8
        )
        assert np.allclose(
            res1.exposures.to_numpy(), res2.exposures.to_numpy(), rtol=1e-6, atol=1e-8
        )

    def test_different_random_state_gives_different_results(self):
        matrix = self.toy_matrix()

        res1 = run_nmf(matrix, n_signatures=2, random_state=0)
        res2 = run_nmf(matrix, n_signatures=2, random_state=1)

        # It's possible (but very unlikely) they end up numerically identical;
        # we test they are *not exactly* equal.
        assert not np.array_equal(
            res1.signatures.to_numpy(), res2.signatures.to_numpy()
        ) or not np.array_equal(
            res1.exposures.to_numpy(), res2.exposures.to_numpy()
        )

    def test_model_params_contain_basic_fields(self):
        matrix = self.toy_matrix()

        res = run_nmf(matrix, n_signatures=2, random_state=0)

        params = res.model_params
        for key in [
            "n_signatures",
            "init",
            "max_iter",
            "random_state",
            "alpha_W",
            "alpha_H",
            "l1_ratio",
            "reconstruction_err_",
            "n_iter_",
        ]:
            assert key in params

class TestSaveLoadNMFResult:
    def _toy_matrix(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "f1": [1.0, 2.0, 3.0],
                "f2": [0.0, 1.0, 0.0],
                "f3": [4.0, 0.0, 1.0],
            },
            index=["s1", "s2", "s3"],
        )

    def test_save_and_load_roundtrip(self, output_dir: str):
        """
        run_nmf -> save_nmf_result -> load_nmf_result
        should round-trip signatures & exposures (up to numeric precision)
        and preserve model_params keys.
        """
        matrix = self._toy_matrix()
        res = run_nmf(matrix, n_signatures=2, random_state=0)

        outdir = Path(output_dir) / "nmf_roundtrip"
        save_nmf_result(res, outdir)

        # Basic files exist
        assert (outdir / "signatures.tsv").is_file()
        assert (outdir / "exposures.tsv").is_file()
        assert (outdir / "metadata.json").is_file()

        res_loaded = load_nmf_result(outdir)

        # Check shapes and labels
        assert list(res_loaded.signatures.index) == list(res.signatures.index)
        assert list(res_loaded.signatures.columns) == list(res.signatures.columns)
        assert list(res_loaded.exposures.index) == list(res.exposures.index)
        assert list(res_loaded.exposures.columns) == list(res.exposures.columns)

        # Numeric close
        assert np.allclose(
            res_loaded.signatures.to_numpy(),
            res.signatures.to_numpy(),
            rtol=1e-6,
            atol=1e-8,
        )
        assert np.allclose(
            res_loaded.exposures.to_numpy(),
            res.exposures.to_numpy(),
            rtol=1e-6,
            atol=1e-8,
        )

        # Metadata contains the core fields
        for key in [
            "n_signatures",
            "init",
            "max_iter",
            "random_state",
            "alpha_W",
            "alpha_H",
            "l1_ratio",
            "format_version",
            "created_at",
        ]:
            assert key in res_loaded.model_params

    def test_load_missing_files_raises(self, output_dir: str):
        """
        If one of signatures.tsv / exposures.tsv / metadata.json is missing,
        load_nmf_result should raise FileNotFoundError.
        """
        outdir = Path(output_dir) / "nmf_missing"
        outdir.mkdir(parents=True, exist_ok=True)

        # Only create one of the files
        (outdir / "signatures.tsv").write_text("dummy\t1\nfeat1\t0.1\n")

        with pytest.raises(FileNotFoundError):
            _ = load_nmf_result(outdir)

class TestProjectOntoSignatures:
    def _toy_matrix(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "f1": [1.0, 2.0, 3.0],
                "f2": [0.0, 1.0, 0.0],
                "f3": [4.0, 0.0, 1.0],
            },
            index=["s1", "s2", "s3"],
        )

    def test_project_basic_nnls_shape_and_nonneg(self):
        """
        Project the same matrix onto its own learned signatures.
        Check:
          - shape (samples x K)
          - non-negativity
          - reconstruction error is smaller than norm of data
        """
        matrix = self._toy_matrix()
        nmf_res = run_nmf(matrix, n_signatures=2, random_state=0)

        exposures_proj = project_onto_signatures(
            new_matrix=matrix,
            signatures=nmf_res.signatures,
            method="nnls",
        )

        # Same samples, same number of signatures
        assert list(exposures_proj.index) == list(matrix.index)
        assert list(exposures_proj.columns) == list(nmf_res.signatures.columns)

        # Non-negative
        arr = exposures_proj.to_numpy()
        assert np.all(arr >= -1e-10)

        # Reconstruction quality: ||X - E_proj * H^T|| < ||X||
        X = matrix.to_numpy(dtype=float)
        H = nmf_res.signatures.to_numpy(dtype=float)  # features x K
        E = exposures_proj.to_numpy(dtype=float)      # samples x K

        X_recon = E @ H.T
        err = np.linalg.norm(X - X_recon, ord="fro")
        base = np.linalg.norm(X, ord="fro")
        assert err < base

    def test_project_with_partial_feature_overlap_warns(self):
        """
        If new_matrix has only a subset of features in signatures.index,
        projection should:
          - warn if overlap < 50%
          - still return a exposures matrix
        """
        matrix = self._toy_matrix()
        nmf_res = run_nmf(matrix, n_signatures=2, random_state=0)

        # new_matrix only has f1
        new_matrix = pd.DataFrame(
            {
                "f1": [1.0, 0.5],
            },
            index=["n1", "n2"],
        )

        with pytest.warns(RuntimeWarning):
            exposures_proj = project_onto_signatures(
                new_matrix=new_matrix,
                signatures=nmf_res.signatures,
                method="nnls",
            )

        assert exposures_proj.shape == (2, nmf_res.signatures.shape[1])
        assert list(exposures_proj.index) == ["n1", "n2"]
        assert list(exposures_proj.columns) == list(nmf_res.signatures.columns)

    def test_project_no_common_features_raises(self):
        """
        If there is no overlap between new_matrix.columns and signatures.index,
        projection should raise ValueError.
        """
        matrix = self._toy_matrix()
        nmf_res = run_nmf(matrix, n_signatures=2, random_state=0)

        new_matrix = pd.DataFrame(
            {
                "g1": [1.0, 2.0],
                "g2": [0.0, 1.0],
            },
            index=["n1", "n2"],
        )

        with pytest.raises(ValueError) as excinfo:
            _ = project_onto_signatures(
                new_matrix=new_matrix,
                signatures=nmf_res.signatures,
                method="nnls",
            )
        assert "No overlapping features" in str(excinfo.value)

    def test_project_invalid_method_raises(self):
        matrix = self._toy_matrix()
        nmf_res = run_nmf(matrix, n_signatures=2, random_state=0)

        with pytest.raises(ValueError):
            _ = project_onto_signatures(
                new_matrix=matrix,
                signatures=nmf_res.signatures,
                method="something-else",
            )
