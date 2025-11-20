# Testing Documentation

## Overview

The `str_mut_signatures` project has a test suite covering:

- VCF parsing and mutation extraction
- Matrix construction and filtering
- NMF decomposition, saving/loading, and projection
- Plotting helpers
- End-to-end pipelines (both Python API and CLI)

**Current Coverage (pytest + coverage):**

| File                                              | Stmts | Miss | Excl | Branches | Partial | Coverage |
|---------------------------------------------------|:-----:|:----:|:----:|:--------:|:-------:|:--------:|
| `src/str_mut_signatures/__init__.py`             |  9    |  0   |  0   |    0     |   0     | 100.00%  |
| `src/str_mut_signatures/cli.py`                  | 166   | 166  |  2   |   54     |   0     |  0.00%   |
| `src/str_mut_signatures/extract_tally/__init__.py` | 0   |  0   |  0   |    0     |   0     | 100.00%  |
| `src/str_mut_signatures/extract_tally/extract_mutations.py` | 117 | 11 | 0 | 48 | 7 | 89.09% |
| `src/str_mut_signatures/extract_tally/filter.py` | 76    |  2   |  0   |   24     |   2     | 96.00%   |
| `src/str_mut_signatures/extract_tally/tally.py`  | 97    |  2   |  0   |   42     |   2     | 97.12%   |
| `src/str_mut_signatures/extract_tally/validate.py` | 79  |  7   |  0   |   28     |   4     | 89.72%   |
| `src/str_mut_signatures/nmf/nmf.py`              | 122   | 13   |  2   |   50     |  15     | 83.72%   |
| `src/str_mut_signatures/nmf/plot.py`             | 100   | 28   |  0   |   42     |  14     | 63.38%   |
| **Total**                                        | 766   | 229  |  4   |  288     |  44     | **68.03%** |

Notes:

- Core extraction, tallying, filtering, and NMF logic are well-covered.
- The CLI module (`cli.py`) is currently tested via the installed console script (subprocess-based CLI integration tests), so direct line coverage is **0%**.
- Plotting functions are partially covered (they create figures but are not deeply asserted).

---

## Quick Start

Run the test suite from the project root:

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=src/str_mut_signatures --cov-report=term-missing

# Run only integration tests (Python API)
pytest tests/integration/ -v

# Run only CLI tests
pytest tests/cli/ -v

# Run unit-style tests in the top-level test files
pytest tests/test_validate.py -v
pytest tests/test_extract_mutations.py -v
pytest tests/test_tally.py -v
pytest tests/test_filter.py -v
pytest tests/test_nmf.py -v
````

Useful switches:

```bash
# Verbose output
pytest tests/ -v

# Stop at first failure
pytest tests/ -x

# Re-run only last failed tests
pytest tests/ --lf

# Show print/log output and locals on failure
pytest tests/ -s -l
```

---

## Test Structure

The test layout is roughly:

```text
tests/
├── conftest.py                     # Global fixtures (data_dir, output_dir, vcf_dir, etc.)
├── data/
│   ├── vcfs/
│   │   └── test_input.zip          # Multiple test VCFs for pipelines
│   ├── pindel_header.vcf           # Header-only VCF for validation tests
│   ├── mutec2_indel.vcf.gz         # VCF without STR annotations (negative test)
│   └── ...                         # Snapshot manifests for integration tests
├── cli/
│   └── test_cli_commands.py        # CLI unit + integration tests
├── integration/
│   └── test_pipeline.py            # Python-API integration + snapshot tests
├── test_validate.py                # Tests for VCF validation logic
├── test_extract_mutations.py       # Tests for parse_vcf_files / process_vcf_to_rows
├── test_tally.py                   # Tests for build_mutation_matrix & helpers
├── test_filter.py                  # Tests for filter_mutation_matrix & FilterSummary
└── test_nmf.py                     # Tests for NMFResult, run_nmf, save/load, projection, plotting
```

This may evolve as the project grows, but the main categories are:

* **Core logic tests**: validation, extraction, tally, filtering, NMF.
* **Integration tests**: full pipeline from VCF to NMF and projection.
* **CLI tests**: invoke the installed `str_mut_signatures` command via `subprocess`.

---

## Key Fixtures

Defined in `tests/conftest.py`:

```python
@pytest.fixture(scope="session")
def data_dir() -> str:
    """Absolute path to the test data directory."""
    return str(Path(__file__).resolve().parent / "data")


@pytest.fixture(scope="session")
def output_dir(tmp_path_factory) -> str:
    """
    Session-wide temporary directory for test outputs.

    Used by tests that compare file hashes or snapshots.
    """
    tmp_dir = tmp_path_factory.mktemp("output")
    return str(tmp_dir)


@pytest.fixture
def temp_output_dir(tmp_path) -> str:
    """
    Per-test temporary output directory.

    Used by tests that just need a scratch directory and
    don't care about reuse across tests.
    """
    d = tmp_path / "out"
    d.mkdir(exist_ok=True)
    return str(d)


@pytest.fixture(scope="session")
def vcf_dir(data_dir: str, output_dir: str) -> str:
    """
    Unpack tests/data/vcfs/test_input.zip into tests/output/vcfs
    and return the directory containing VCF files.

    If the directory already exists, it is reused.
    """
    ...
```

These fixtures give all tests consistent access to:

* Stable **input data** (`data_dir`).
* A **session-wide output directory** (`output_dir`) for snapshot tests.
* A **scratch per-test directory** (`temp_output_dir`).
* A directory of **test VCFs** (`vcf_dir`) used by extraction and pipeline tests.

---

## Test Categories

### 1. Validation (`test_validate.py`)

Covers `str_mut_signatures.extract_tally.validate`:

* `validate_vcf`:

  * Detects presence/absence of STR annotations (`RU`, `REF`, `REPCN`).
  * Validates that there are at least two samples.
  * Captures normal/tumor sample names and phasing information.

* Tests use small VCFs:

  * `pindel_header.vcf` – valid structure + annotations.
  * `mutec2_indel.vcf.gz` – negative test: no STR annotations.

**Example:**

```bash
pytest tests/test_validate.py -v
```

---

### 2. Mutation Extraction (`test_extract_mutations.py`)

Covers `extract_tally/extract_mutations.py`:

* `process_vcf_to_rows(path)`:

  * Parses one VCF at a time (supports `.vcf` and `.vcf.gz`).
  * Applies optional filtering:

    * `FILTER == "PASS"`.
    * `INFO.PERFECT == TRUE`.
  * Extracts allele-level repeat copy numbers for normal/tumor:

    * `normal_allele_a`, `normal_allele_b`
    * `tumor_allele_a`,  `tumor_allele_b`
  * Captures STR metadata (`end`, `period`, `ref`, `motif`).
  * Stores genotype separator (`genotype_separator`) for phasing.

* `parse_vcf_files(input_dir)`:

  * Walks a directory of VCFs.
  * Aggregates all rows into a single `DataFrame`.
  * Skips files that raise errors, printing helpful messages.
  * Ensures the expected columns are present and that there is at least one row.

**Example:**

```bash
pytest tests/test_extract_mutations.py -v
```

---

### 3. Tally / Matrix Construction (`test_tally.py`)

Covers `extract_tally/tally.py`:

* `validate_mutations_data(df)`:

  * Ensures required columns are present.
  * Resolves motif column (`motif` vs `RU`).
  * Detects whether `genotype_separator` is available.

* `compute_changes_for_row(row)`:

  * Phased (GT with `|`): computes allele-specific deltas and refs.
  * Unphased (`/` or missing): computes a single combined event.

* `motif_is_at_rich(motif)`:

  * Classifies motifs as `AT_rich` vs `non_AT_rich`.

* `make_feature(...)`:

  * Builds feature keys based on:

    * `ru` mode (`None`, `"length"`, `"ru"`, `"AT"`),
    * `ref_length`,
    * `change`.

* `build_mutation_matrix(...)`:

  * Converts long-format mutation data to a count matrix (samples × features).
  * Properly handles phased/unphased genotypes.
  * Supports all `ru` modes, `ref_length`, and `change` combinations.
  * Drops non-somatic events if `change=True`.

**Example:**

```bash
pytest tests/test_tally.py -v
```

---

### 4. Filtering (`test_filter.py`)

Covers `extract_tally/filter.py`:

* `filter_mutation_matrix(matrix, feature_method=..., ...)`:

  * `feature_method="manual"`:

    * Uses `min_feature_total`, `min_samples_with_feature`, `min_sample_total`.

  * `feature_method="elbow"`:

    * Finds an elbow in sorted feature totals.
    * Uses the inferred threshold + the same sample constraints.

  * `feature_method="percentile"`:

    * Keeps features above a configurable percentile of feature totals.

* Returns:

  * Filtered matrix.
  * `FilterSummary` with counts of samples/features before/after filtering and thresholds used.

Tests also exercise:

* Edge cases (small matrices, all-zero rows/columns, tight thresholds).
* The safety of the elbow calculation (no division by zero, no crashes).

**Example:**

```bash
pytest tests/test_filter.py -v
```

---

### 5. NMF & Projection (`test_nmf.py`)

Covers `nmf/nmf.py` and `nmf/plot.py`:

* `NMFResult` dataclass structure.

* `run_nmf`:

  * Validates non-negativity.
  * Preserves indices/columns.
  * Stores hyperparameters and diagnostics in `model_params`.

* `save_nmf_result` / `load_nmf_result`:

  * Round-trip equality of shapes and labels.
  * Metadata consistency (`format_version`, shapes, column names, etc.).

* `project_onto_signatures`:

  * Handles feature intersection between new matrices and learned signatures.
  * Uses NNLS (`scipy.optimize.nnls`) when available.
  * Ensures exposures are non-negative and shaped correctly.

* Plotting:

  * `plot_signatures(nmf_res, top_n=...)` – returns a `matplotlib` Figure.
  * `plot_exposures(nmf_res)` – returns a Figure.
  * `compute_pca` and `plot_pca_samples` – tested for shape and basic behaviour.

**Example:**

```bash
pytest tests/test_nmf.py -v
```

---

### 6. Integration Pipeline (`tests/integration/test_pipeline.py`)

End-to-end test using the **Python API**:

1. `parse_vcf_files(vcf_dir)` – load multiple VCFs.
2. `build_mutation_matrix(...)` – construct a matrix.
3. `filter_mutation_matrix(...)` – run all filter methods.
4. `run_nmf(...)` – fit an NMF model.
5. Plot signatures, exposures, PCA (just asserting no crash).
6. `save_nmf_result` and `load_nmf_result` – round-trip.
7. `project_onto_signatures` on a new single-VCF matrix.
8. Build a **hash manifest** of all files produced in a specific output subdirectory (e.g. `integration_pipeline/`) and compare against a golden snapshot (`tests/data/...txt`).

Snapshot logic:

* All files under a pipeline output directory are hashed via SHA-256.

* Paths and hashes are written as:

  ```text
  relative/path/to/file\t<sha256>
  ```

* First run can bootstrap the golden file; subsequent runs ensure outputs are stable.

**Example:**

```bash
pytest tests/integration/test_pipeline.py -v
```

---

### 7. CLI Integration (`tests/cli/test_cli_commands.py`)

End-to-end test using the **installed CLI** (`str_mut_signatures`):

* Invokes subcommands via `subprocess.run`:

  * `str_mut_signatures extract ...`
  * `str_mut_signatures filter ...`
  * `str_mut_signatures nmf ...`
  * `str_mut_signatures project ...`

* Asserts:

  * Exit code 0 on success.
  * Meaningful error codes and messages for invalid paths/options.
  * Required output files exist (`matrix_*.tsv`, `signatures.tsv`, `exposures.tsv`, `metadata.json`, new exposures).

* Builds a hash manifest (similar to the Python-API integration test) for a CLI-specific output directory (e.g. `cli_integration_pipeline/`) and compares it against a golden snapshot (`tests/data/test_cli_full_pipeline_from_vcf_to_projection_and_snapshot.txt`).

**Example:**

```bash
pytest tests/cli/test_cli_commands.py -v
```

---

## Running With Coverage

```bash
# Full coverage report in terminal (missing lines)
pytest tests/ --cov=src/str_mut_signatures --cov-report=term-missing

# HTML coverage (open htmlcov/index.html afterwards)
pytest tests/ --cov=src/str_mut_signatures --cov-report=html
```

Notes:

* `cli.py` shows 0% coverage because CLI tests run the installed console script rather than importing `str_mut_signatures.cli` directly. If you want direct line coverage there, add tests that import and call `cli.main()`.

---

## Debugging Tests

Handy pytest options:

```bash
# Show print/log output
pytest tests/ -s

# Stop at first failure
pytest tests/ -x

# Show local variables on failure
pytest tests/ -l

# Only re-run previously failing tests
pytest tests/ --lf

# Drop into pdb on failure
pytest tests/ --pdb
```

---

## Guidelines for New Tests

1. **Use descriptive names**

   * ✅ `test_filters_low_count_features()`
   * ✅ `test_run_nmf_raises_on_negative_values()`
   * ❌ `test_function1()`

2. **Test one concept per test**
   Keep each test focused (e.g. one aspect of filtering, one error condition).

3. **Use fixtures for setup**
   Reuse `data_dir`, `output_dir`, `vcf_dir`, and `temp_output_dir` instead of duplicating setup code.

4. **Test edge cases**

   * Empty matrices
   * Single-sample / single-feature matrices
   * VCFs missing required annotations
   * Extremely sparse matrices for NMF/projection

5. **Test error conditions**

   * Missing or invalid paths in CLI.
   * Negative entries for NMF.
   * Projection with no overlapping features.

6. **Keep integration tests stable**

   * If you intentionally change the output format, update the golden hash manifests after verifying that the changes are expected:

     * Python API: `tests/data/test_full_pipeline_from_vcf_to_projection_and_snapshot.txt`
     * CLI: `tests/data/test_cli_full_pipeline_from_vcf_to_projection_and_snapshot.txt`

---

**Last Updated**: based on test suite and coverage as of the current `str_mut_signatures` development state (overall coverage ~68%).
