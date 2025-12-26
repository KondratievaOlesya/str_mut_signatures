# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
as far as reasonably possible for a research codebase.
## [1.0.0] - 2025-12-22
### Added

* **Persistent sample grouping in `NMFResult`**

  * Introduced a new `groups` attribute in `NMFResult`:

    * Stored as a `pandas.DataFrame` with:

      * index = sample IDs
      * column = `"group"`
    * Always present (at minimum, all samples assigned to group `"1"`).

* **Optional clustering at NMF runtime**

  * `run_nmf()` now supports optional sample clustering via `max_clusters`:

    * If `max_clusters > 1`, samples are clustered based on **signature exposures** using KMeans.
    * Optimal number of clusters is selected using **silhouette score**.
    * Cluster assignments are stored in `NMFResult.groups`.
    * If clustering is disabled or not possible, all samples are assigned to group `"1"`.

* **Group-aware PCA plotting**

  * `plot_pca_samples()` now:

    * Colors samples **exclusively using `NMFResult.groups`**
    * Automatically aligns PCA coordinates and group labels by sample ID
     * Automatically detects whether `group` is **continuous** vs **categorical**:
        * Treats `group` as continuous only if it is numeric **and** has sufficient diversity
          (more than **10** unique non-null values **or** > **30%** unique fraction)
        * Continuous groups are shown with a **color gradient + colorbar**
        * Categorical groups are shown with **separate legend entries**
    * Does **not** perform clustering internally
    * Returns:

      * PCA coordinates
      * explained variance ratio
      * matplotlib axes

* **Group-aware exposure plotting**

  * `plot_exposures()` now:

    * Uses `NMFResult.groups["group"]` for:

      * sample ordering
      * visual separation between groups
    * No longer performs clustering internally
    * Ensures consistent grouping across all downstream plots

* **Groups persisted in I/O**

  * `save_nmf_result()` now writes:

    * `groups.tsv`
  * `load_nmf_result()` restores `groups` alongside signatures and exposures.

### Changed

* **Plotting API semantics**

  * Clustering logic has been **fully decoupled from plotting**:

    * Clustering is now a modeling decision (`run_nmf`)
    * Plotting functions are purely representational
  * Removed `cluster` / `color_by` logic from PCA and exposure plotting.

* **End-to-end pipeline consistency**

  * Integration tests updated to:

    * Validate presence and correctness of `groups`
    * Use group-aware plotting APIs
    * Ensure consistent sample alignment across NMF, PCA, and exposure plots

### Removed

* **Internal clustering from plotting functions**

  * `plot_pca_samples()` and `plot_exposures()` no longer:

    * run KMeans
    * select `k`
    * return cluster labels
  * All grouping information must come from `NMFResult.groups`.

## [0.3.0] - 2025-12-02

### Added

- **Motif normalization**

  - Introduced `normalize_motif()` to standardize all repeat-unit (RU) motifs:

    - Converts motifs to **uppercase**.
    - Removes spaces and any non-letter characters.
    - Canonicalizes motifs using **`GetCanonicalMotif`** from `trtools.utils.utils`.

- **Support for multiple STR callers (GangSTR and conSTRain)**  
    The VCF validator and mutation extractor now recognize and parse STR annotations from both tools:
    - GangSTR (`FORMAT/REPCN`)
    - conSTRain (`FORMAT/REPLEN`)
        
    
    The parser automatically detects:
    - caller type (`gangstr`, `constrain`, or `unknown`)
    - correct copy-number field (`REPCN` or `REPLEN`)
    - STR annotation presence via `INFO/RU` and `INFO/REF`.
        
- **Automatic extraction of STR allele counts**  
    A unified `parse_copy_number()` replaces the old REPCN-specifi `parse_repcn` parser.  
        
- **New integration test** for real conSTRain VCF
        
- **Clustering support for PCA and exposure plots**  
    A unified mechanism:
    - automatically selects optimum number of clusters (`k`) using silhouette score
    - groups samples visually and orders them within clusters
    - exposes cluster labels in results
        
- **Enhanced exposure plotting (`plot_exposures`)**
    - Pagination for large sample sets (`max_samples_per_fig`)
    - Cluster-aware ordering (with visible gaps between clusters)
    - Dual-plot mode:
        - **absolute exposures**
        - **proportional exposures** (stacked, normalized to 1)
            
- **Enhanced PCA plotting API**
    - `plot_pca_samples()` is now the main high-level function:
        - accepts `NMFResult` directly
        - computes PCA internally
        - optional clustering with cluster-based coloring
            
    - Returns:
        - PCA coordinates
        - variance explained
        - cluster labels
        - matplotlib axes
            

### **Deprecated**
- `parse_repcn()`  
    Replaced by the more general `parse_copy_number()`.

## [0.2.1] - 2025-11-20

- Dependencies package fix

## [0.2.0] - 2025-11-20

### Added

- **New extraction / tally architecture**

  - Introduced `str_mut_signatures.extract_tally` subpackage with:
    - `extract_mutations.py` for streaming per-VCF mutation extraction.
    - `tally.py` for building mutation count matrices.
    - `validate.py` for lightweight VCF validation.

  - New `VCFValidationResult` (`NamedTuple`) returned by `validate_vcf`, summarising:
    - Presence of STR annotations (`RU`, `REF`, `REPCN`).
    - Whether there are at least two samples.
    - Normal/tumor sample names.
    - Whether genotypes are phased or unphased (via GT separator).

- **VCF parsing + phasing support**

  - `parse_vcf_files(vcf_dir)` now:
    - Processes VCF files one at a time (memory-efficient).
    - Supports both plain and gzipped VCFs.
    - Optionally filters by `FILTER == PASS` and `INFO.PERFECT == TRUE`.
    - Extracts per-sample allele repeat counts (`REPCN`) into:
      - `normal_allele_a`, `normal_allele_b`
      - `tumor_allele_a`,  `tumor_allele_b`
    - Tracks genotype phasing via a `genotype_separator` column (`'|'` or `'/'`).

- **More expressive mutation matrix construction**

  - `build_mutation_matrix` has been rewritten around small helpers:
    - `validate_mutations_data` – checks required columns and motif column name.
    - `compute_changes_for_row` – computes allele-level or combined changes:
      - Phased (`'|'`): two allele-specific changes per locus.
      - Unphased or missing (`'/'`/None): a single combined change per locus.
    - `motif_is_at_rich` – classifies motifs as `AT_rich` vs `non_AT_rich`.
    - `make_feature` – constructs feature keys consistently.

  - New `RuMode` type and new `ru` mode:
    - `ru=None` – no motif information.
    - `ru="length"` – motif length: `LEN1`, `LEN2`, …
    - `ru="ru"` – full motif string: `A`, `AT`, `AAT`, …
    - **`ru="AT"` – new**: motif-level classification:
      - `AT_rich`    (only A/T).
      - `non_AT_rich` (contains C/G).

  - Matrix logic now:
    - Correctly handles phased vs unphased genotypes.
    - Uses normal alleles as reference length proxy (per-allele or combined).
    - Optionally includes `ref_length` and `change` in the feature key.
    - Operates in a long format internally (one event per sample/allele),
      then pivots back to a samples × features count matrix.

- **Filtering of sparse mutation matrices**

  - New `filter_mutation_matrix` in `str_mut_signatures.extract_tally.filter`:
    - Supports multiple strategies:

      - `feature_method="manual"`:
        - Keep features with total counts ≥ `min_feature_total`.
        - Require presence in ≥ `min_samples_with_feature` samples.
        - Drop samples with totals < `min_sample_total`.

      - `feature_method="elbow"`:
        - Use an “elbow” heuristic on feature totals to choose a cutoff.
        - Apply the same minimum samples and sample total constraints.

      - `feature_method="percentile"`:
        - Keep features above a specified percentile of feature totals
          (e.g. `feature_percentile=0.9`).

    - Returns:
      - Filtered matrix.
      - A `FilterSummary` dataclass with basic stats (e.g. counts before/after,
        thresholds used).

- **NMF module with reusable model container**

  - New `str_mut_signatures.nmf.nmf` module using `sklearn.decomposition.NMF`.

  - Added `NMFResult` dataclass:
    - `signatures`: `DataFrame` (features × K).
    - `exposures`: `DataFrame` (samples × K).
    - `model_params`: dict of hyperparameters and metadata
      (`n_signatures`, `init`, `max_iter`, `random_state`, etc.).

  - New `run_nmf` function:
    - Validates that the input matrix is non-negative.
    - Preserves row and column labels through the NMF.
    - Exposes key parameters:
      - `n_signatures`, `init`, `max_iter`, `random_state`,
      - `alpha_W`, `alpha_H`, `l1_ratio`.
    - Stores reconstruction error, number of iterations, and a simple
      `format_version` in `model_params`.

- **Saving and loading NMF results**

  - `save_nmf_result(result, outdir)`:
    - Writes:
      - `signatures.tsv`   (features × K).
      - `exposures.tsv`    (samples × K).
      - `metadata.json`    (hyperparameters, shapes, column names,
                            `format_version`, `created_at` timestamp, etc.).

  - `load_nmf_result(outdir) -> NMFResult`:
    - Reconstructs `NMFResult` from TSVs and JSON.
    - Ensures signatures/exposures indices and columns are restored.

- **Projecting new samples onto existing signatures**

  - Added `project_onto_signatures(new_matrix, signatures, method="nnls")`:
    - Aligns features between `new_matrix` and `signatures` via intersection.
    - Uses non-negative least squares (`scipy.optimize.nnls`) per sample
      to compute exposures for existing signatures.
    - Returns a new exposures matrix (samples × K).
    - Gracefully degrades / can be skipped if `scipy` is not installed.

- **Plotting utilities**

  - New `str_mut_signatures.nmf.plot` module:

    - `plot_signatures(NMFResult, top_n=20, ...)`:
      - Simple barplots for each signature, highlighting top contributing
        features.

    - `plot_exposures(NMFResult, sort_by=None, ...)`:
      - Stacked or grouped barplots of signature exposures per sample.

    - `compute_pca(matrix_or_exposures, n_components=2)`:
      - Returns `(coords_df, explained_variance_ratio)`, where `coords_df`
        has columns `PC1`, `PC2`.

    - `plot_pca_samples(coords, labels=None, ...)`:
      - Scatter plot of samples in PC space (e.g. exposures PCA).

- **Richer CLI with full pipeline support**

  - CLI now has four subcommands:

    - `extract`:
      - Reads a directory of VCFs and writes a count matrix.

    - **`filter` (new)**:
      - Runs `filter_mutation_matrix` with manual, elbow, or percentile
        strategies.
      - Writes a filtered matrix.

    - `nmf`:
      - Runs NMF on a matrix and saves results using `save_nmf_result`
        (signatures, exposures, metadata).

    - **`project` (new)**:
      - Loads an existing NMF result directory.
      - Projects a new count matrix onto learned signatures and writes
        new exposures.

  - CLI options now expose:
    - `ru` mode including `"AT"`.
    - Filter thresholds/percentiles.
    - NMF parameters (`n_signatures`, `max_iter`, `random_state`, `alpha_W`,
      `alpha_H`, `l1_ratio`).

- **Integration tests & output snapshotting**

  - Added end-to-end integration tests for:

    - Library API pipeline:
      - VCF → parse → matrix → filter → NMF → save/load → project.

    - CLI pipeline:
      - `extract → filter → nmf → project` via the actual `str_mut_signatures`
        console script.

  - Tests now:
    - Save all intermediate outputs into `tests/output/...`.
    - Build a **hash manifest** (relative path + SHA-256) over output dirs.
    - Compare against golden manifest files in `tests/data/` to detect
      any regression in generated outputs.
    - Ignore `metadata.json` in manifests to avoid timestamp-based noise.

### Changed

- **Top-level NMF API**

  - `run_nmf` (imported from `str_mut_signatures`) now returns an `NMFResult`
    object instead of a plain `(W, H)` tuple.
  - Code that previously expected `(W, H)` should be updated to:

    ```python
    nmf_res = run_nmf(matrix, n_signatures=K)
    signatures = nmf_res.signatures
    exposures = nmf_res.exposures
    ```

- **CLI NMF behaviour**

  - The `nmf` subcommand now:
    - Uses the new `run_nmf` + `save_nmf_result` pipeline.
    - Produces `signatures.tsv`, `exposures.tsv`, and `metadata.json`
      in the output directory.
  - Logging messages and argument validation have been updated accordingly.

- **Mutation matrix computation**

  - `build_mutation_matrix`:
    - Now explicitly supports phased vs unphased genotypes via
      `genotype_separator`.
    - Uses per-allele changes when phased; otherwise uses a single combined
      event per locus.
    - Builds feature keys through the new `make_feature` helper and supports
      the new `"AT"` motif mode.

- **CLI error messages**

  - Validation errors for missing files/dirs now use consistent “not found”
    wording, e.g.:
    - `"--matrix file not found: /path/to/file.tsv"`.
  - Errors are logged as `ERROR - Validation error: ...` to satisfy existing
    tests that look for “not found” / “no such file” patterns.

- **Filter CLI logging**

  - The CLI no longer relies on specific fields such as `n_features_kept`
    on the `FilterSummary` object.
  - Instead, it reports shapes based on the input and filtered matrices
    directly (robust to future changes in the summary dataclass).

### Fixed

- More robust VCF parsing:

  - Malformed lines and entries with non-numeric `REPCN` values are safely
    skipped.
  - Gzipped and plain VCFs are handled transparently.
  - Normal/tumor sample ordering is validated; meaningful errors are raised
    if assumptions are violated.

- Integration tests:

  - Snapshots ignore `metadata.json` to avoid failures caused by
    date/time metadata changes.
  - CLI integration now handles filter behaviour correctly after
    decoupling from `FilterSummary` internals.

---

## [0.1.0] - 2025-11-04

_Initial release_

### Added

- Basic support for:

  - Parsing STR-annotated paired tumor–normal VCFs.
  - Constructing mutation count matrices using motif length, reference length,
    and tumor–normal change.
  - A simple NMF wrapper and minimal CLI (`extract`, `nmf`).
