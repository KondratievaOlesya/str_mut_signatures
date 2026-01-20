==================
str_mut_signatures
==================

.. image:: https://img.shields.io/badge/docs-GitHub%20Pages-blue
   :target: https://acg-team.github.io/str_mut_signatures/

.. image:: https://github.com/acg-team/str_mut_signatures/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/acg-team/str_mut_signatures/actions/workflows/ci.yml
   :alt: Tests

.. image:: https://codecov.io/gh/acg-team/str_mut_signatures/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/acg-team/str_mut_signatures
   :alt: Coverage

STR Mutation Signature Analysis.

Python package for analysis of Short Tandem Repeat (STR) mutation signatures from VCF files.
It extracts somatic STR mutation events from paired tumor–normal VCFs, builds
count matrices, filters them, performs NMF-based signature decomposition, and
projects new samples onto learned STR mutation signatures.

Contents
========

- `Installation`_
- `Quick start`_
- `Input format`_
- `Annotating standard VCFs`_
- `Matrix construction`_
- `Filtering mutation matrices`_
- `NMF signatures and projection`_
- `Visualization`_
- `Command line interface`_
- `Python API`_
- `Output`_
- `Contributing`_
- `License`_


Installation
============

From PyPI
---------

The package is available through `PyPI <https://pypi.org/project/str_mut_signatures/>`_. Install with:

.. code-block:: shell

    pip install str_mut_signatures


From source
-----------

.. code-block:: shell

    git clone https://github.com/acg-team/str_mut_signatures
    cd str_mut_signatures
    pip install -e .

Development installation
------------------------

.. code-block:: shell

    pip install -r requirements_dev.txt


Quick start
===========

Python Library Usage
--------------------

Basic pipeline
~~~~~~~~~~~~~~

.. code-block:: python

    from str_mut_signatures import (
        parse_vcf_files,
        build_mutation_matrix,
        filter_mutation_matrix,
        run_nmf,
        save_nmf_result,
        load_nmf_result,
        project_onto_signatures,
        plot_exposures,
        plot_pca_samples,
        plot_signatures,
    )

    # 1) Parse annotated paired tumor–normal VCF files into a long table
    mutations = parse_vcf_files("vcf_directory/")

    # 2) Build a mutation count matrix
    # ru_length:
    #   include repeat-unit length (LEN{len(motif)})
    # ru:
    #   None     -> ignore motif
    #   "ru"     -> use full repeat unit sequence (e.g. AT, AAT)
    #   "class"  -> AT-only, GC-only or mixed classification
    # ref_length:
    #   include reference repeat length as a feature component
    # change:
    #   include tumor–normal repeat-length change
    matrix = build_mutation_matrix(
        mutations,
        ru_length=True,
        ru="class",
        ref_length=True,
        change=True,
    )

    # 3) Filter the matrix (e.g. manual thresholds)
    matrix_filt, summary = filter_mutation_matrix(
        matrix,
        feature_method="manual",
        min_feature_total=10,
        min_samples_with_feature=3,
        min_sample_total=0,
    )

    # 4) Run NMF
    # Optional clustering is performed here and stored in nmf_res.groups
    nmf_res = run_nmf(
        matrix_filt, 
        n_signatures=5, 
        random_state=11, 
        max_clusters=4 # <= 1 disables clustering
    )

    # Access signatures and exposures
    signatures = nmf_res.signatures   # features x K
    exposures  = nmf_res.exposures    # samples x K
    groups     = nmf_res.groups       # samples x 1 ("group")

    # 4a) Plot mutation signatures
    fig = plot_signatures(nmf_res, top_n=20)

    # 4b) Visualize sample exposures to signatures (ordered by group, then total exposure)
    figs = plot_exposures(
        nmf_res,
        stacked=True,
        max_samples_per_fig=40,  # paginate if many samples
    )

    # figs["absolute"]   -> list of figures with absolute exposures
    # figs["proportion"] -> list of figures with exposures normalized per sample

    # 4c) Visualize samples in PCA space based on exposures (colored by stored groups)
    coords, var_ratio, ax = plot_pca_samples(nmf_res)

    # coords   -> DataFrame with PC1, PC2, ...
    # var_ratio-> variance explained by PC1 and PC2
    # ax       -> matplotlib Axes with the scatter plot

    # 5) Save NMF result for reuse
    save_nmf_result(nmf_res, "nmf_results")

    # 6) Load NMF result and project new samples
    nmf_loaded = load_nmf_result("nmf_results")
    new_exposures = project_onto_signatures(
        new_matrix=new_counts_df,
        signatures=nmf_loaded.signatures,
        method="nnls",
    )

Command Line
------------

1. **Extract** somatic STR mutation counts from paired tumor–normal VCFs:

.. code-block:: shell

    str_mut_signatures extract \
        --vcf-dir data/vcfs/ \
        --out-matrix counts_raw.tsv \
        --ru-length \
        --ru class \
        --ref-length \
        --change

This produces a count matrix (TSV) with:

- rows = samples
- columns = STR mutation features such as:

.. code-block:: text

    LEN{motif_length}_{motif class}_{ref_length}_{change}

For example: ``LEN1_AT_only_10_+1`` means:

- motif length = 1 bp
- motif base class = AT-only
- reference repeat length = 10 copies
- tumor has +1 copy relative to normal.

2. **Filter** the matrix to remove extremely rare features/samples:

.. code-block:: shell

    str_mut_signatures filter \
        --matrix counts_raw.tsv \
        --out-matrix counts_filtered.tsv \
        --feature-method elbow

3. **Run NMF** to learn STR mutation signatures:

.. code-block:: shell

    str_mut_signatures nmf \
        --matrix counts_filtered.tsv \
        --outdir nmf_results \
        --n-signatures 5

This writes:

- ``nmf_results/signatures.tsv`` – STR mutation signatures (features x K)
- ``nmf_results/exposures.tsv`` – sample exposures (samples x K)
- ``nmf_results/metadata.json`` – parameters and metadata.

4. **Project new samples** onto existing signatures:

.. code-block:: shell

    str_mut_signatures project \
        --matrix new_counts.tsv \
        --nmf-dir nmf_results \
        --out-exposures new_exposures.tsv

Input format
============

To be processed by ``str_mut_signatures``, VCF files must:

1. Contain **paired samples** (normal and tumor) per record.
2. Be annotated with STR-specific fields that describe the repeat motif and
   **allele-level repeat copy numbers**.

Supported sources
-----------------

``str_mut_signatures`` currently supports STR VCFs produced by:

- `GangSTR  <https://github.com/gymreklab/GangSTR>`_
- `conSTRain  <https://github.com/acg-team/ConSTRain>`_
- VCFs annotated with the `strvcf_annotator  <https://github.com/acg-team/strvcf_annotator/>`_ tool

as long as they satisfy the requirements below.

Required structure
------------------

- Each VCF record must contain at least two samples:

  - **Sample 1 (first column after FORMAT): normal**
  - **Sample 2 (second column after FORMAT): tumor**

  By default, ``str_mut_signatures`` assumes this order and computes **somatic**
  changes as tumor vs normal. Only loci with differences between tumor and normal
  are used (somatic STR mutations).

- Required annotations:

  **INFO fields**

  - ``RU``: Repeat unit / motif (e.g. ``A``, ``AT``, ``AAT``).
  - ``REF``: Reference repeat count (copy number of the motif in the reference genome).

  **FORMAT fields (copy number)**

  One of the following per sample:

  - ``REPCN``: Per-allele repeat copy number (used by GangSTR and ``strvcf_annotator``).
  - ``REPLEN``: Per-allele repeat length in repeat units (used by conSTRain).

  ``str_mut_signatures`` automatically detects which field is present and uses it
  as the source of allele copy numbers.

Example schemas
---------------

**GangSTR / strvcf_annotator-style**

.. code-block:: text

    ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
    ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
    ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
    #CHROM POS  ID REF ALT QUAL FILTER INFO        FORMAT        NORMAL     TUMOR
    chr1   100 .  A   AT  .    .     RU=A;REF=10  GT:REPCN      0/0:10,10  0/1:10,11

**conSTRain-style**

.. code-block:: text

    ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
    ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
    ##FORMAT=<ID=REPLEN,Number=R,Type=Integer,Description="Per-allele repeat length in repeat units">
    #CHROM POS  ID REF ALT QUAL FILTER INFO        FORMAT        NORMAL       TUMOR
    chr1   100 .  A   AT  .    .     RU=A;REF=10  GT:REPLEN     0/0:10,10    0/1:10,11

From such records, ``str_mut_signatures``:

- Extracts **allele-level repeat copy numbers** from ``REPCN`` or ``REPLEN``.
- Compares NORMAL vs TUMOR allele counts at each locus.
- Identifies loci where tumor repeat copy number differs from normal.
- Encodes the net repeat-length **change** as tumor–normal (e.g. ``+1``).
- Uses only these somatic STR events for downstream count matrices and signature extraction.


Annotating standard VCFs
========================

If your VCFs lack ``RU``, ``REF``, or copy-number fields, you can annotate them using
the tool ``strvcf_annotator``:

- Takes standard VCF + STR reference.
- Produces STR-annotated VCFs compatible with ``str_mut_signatures``.

For details see: `strvcf_annotator  <https://github.com/acg-team/strvcf_annotator/>`_.


Matrix construction
===================

``build_mutation_matrix`` defines feature keys using independent options.

Core components
---------------

Each feature key can include:

- repeat-unit length (``LEN{len(motif)}``) via ``ru_length``
- repeat-unit content via ``ru`` (optional)
- reference repeat count via ``ref_length``
- tumor–normal delta via ``change``

Repeat-unit content modes
-------------------------

- ``ru`` omitted (``ru=None``):
  no repeat-unit content is included.

- ``ru="class"``:
  base composition class is used:

  - ``AT_only``: motif contains only A/T
  - ``GC_only``: motif contains only G/C
  - ``mixed``: mixed A/T and G/C

- ``ru="ru"``:
  full repeat-unit sequence is used (e.g. ``A``, ``AT``, ``AAT``).

Examples
--------

1. Motif length + base class + ref length + somatic change:

.. code-block:: python

    m = build_mutation_matrix(
        mutations,
        ru_length=True,
        ru="class",
        ref_length=True,
        change=True,
    )
    # e.g. LEN2_AT_only_15_+1

2. Full motif + ref length + somatic change:

.. code-block:: python

    m = build_mutation_matrix(
        mutations,
        ru_length=False,
        ru="ru",
        ref_length=True,
        change=True,
    )
    # e.g. AT_15_+1

3. Length only:

.. code-block:: python

    m = build_mutation_matrix(
        mutations,
        ru_length=True,
        ru=None,
        ref_length=False,
        change=True,
    )
    # e.g. LEN1_+1

4. Presence/absence-style summaries (no delta term):

.. code-block:: python

    m = build_mutation_matrix(
        mutations,
        ru_length=True,
        ru=None,
        ref_length=True,
        change=False,
    )
    # e.g. LEN1_10




Filtering mutation matrices
===========================

Large STR feature spaces can be sparse. ``filter_mutation_matrix`` provides
several strategies to reduce noise before NMF.

Supported methods
-----------------

.. code-block:: python

    from str_mut_signatures import filter_mutation_matrix

    filtered, summary = filter_mutation_matrix(
        matrix,
        feature_method="manual",
        min_feature_total=10,
        min_samples_with_feature=3,
        min_sample_total=0,
        feature_percentile=0.9,  # used for percentile method
    )

Methods:

- ``feature_method="manual"``

  - Keep features with:

    - total count across samples >= ``min_feature_total``
    - present (non-zero) in at least ``min_samples_with_feature`` samples.

  - Drop samples with total counts < ``min_sample_total``.

- ``feature_method="elbow"``

  - Compute feature totals.
  - Use an "elbow" heuristic to choose a count threshold.
  - Keep features above that threshold.
  - Apply ``min_samples_with_feature`` and ``min_sample_total`` as in manual mode.

- ``feature_method="percentile"``

  - Compute feature totals.
  - Keep features above a chosen percentile of totals
    (e.g. ``feature_percentile=0.9`` keeps the top 10% by total count).
  - Apply ``min_samples_with_feature`` and ``min_sample_total`` as in manual mode.

The function returns:

- ``filtered`` – filtered count matrix.
- ``summary`` – small dataclass with filtering statistics (e.g. numbers of features/samples before/after, thresholds used).


NMF signatures and projection
=============================

NMF decomposition
-----------------

NMF is used to decompose the filtered matrix into:

- **Signatures**: STR mutation patterns (features x K)
- **Exposures**: how much each sample uses each signature (samples x K)

.. code-block:: python

    from str_mut_signatures import run_nmf

    nmf_res = run_nmf(
        matrix,
        n_signatures=5,
        init="nndsvd",
        max_iter=200,
        random_state=0,
        alpha_W=0.0,
        alpha_H=0.0,
        l1_ratio=0.0,
    )

    signatures = nmf_res.signatures  # DataFrame: features x K
    exposures  = nmf_res.exposures   # DataFrame: samples x K
    params     = nmf_res.model_params

Saving and loading NMF results
------------------------------

You can save and reload NMF results in a stable format (TSV + JSON):

.. code-block:: python

    from str_mut_signatures import save_nmf_result, load_nmf_result

    save_nmf_result(nmf_res, "nmf_results")

    nmf_loaded = load_nmf_result("nmf_results")
    # nmf_loaded.signatures, nmf_loaded.exposures, nmf_loaded.model_params

Projecting new samples
----------------------

Given a previously learned set of signatures, you can compute exposures for
new samples (e.g. a new cohort or single sample):

.. code-block:: python

    from str_mut_signatures import project_onto_signatures

    new_exposures = project_onto_signatures(
        new_matrix=new_counts_df,
        signatures=nmf_loaded.signatures,
        method="nnls",  # non-negative least squares
    )

Rows in ``new_exposures`` are new samples, columns are signatures.


Visualization
=============

Overview
--------

After fitting NMF, ``str_mut_signatures`` provides convenience functions to:

- **Visualize signature exposures per sample** (bar plots).
- **Visualize samples in PCA space** (scatter plots of PC1 vs PC2).

These plots are useful for:

- Checking whether signatures separate known biological groups.
- Detecting outlier samples or batch effects.
- Inspecting clusters of samples with similar STR mutation patterns.
- Communicating results in figures for manuscripts or presentations.

Exposure bar plots
------------------

Use ``plot_exposures`` to visualize how much each sample uses each signature.

Key features:

- Stacked or grouped bars.
- Optional pagination when there are many samples.
- Samples are ordered by:
  1. group (ascending)
  2. total exposure within each group (descending)

.. code-block:: python

    from str_mut_signatures import plot_exposures

    # nmf_res is the result of run_nmf(...)
    figs = plot_exposures(
        nmf_res,
        stacked=True,
        max_samples_per_fig=40,  # paginate if many samples
    )

    # figs["absolute"]   -> list of figures with absolute exposures
    # figs["proportion"] -> list of figures with exposures normalized per sample

Typical use:

- **Absolute exposures**: how much each signature contributes in raw counts.
- **Proportions**: composition of each sample (bars sum to 1), easier for comparing
  relative contributions across samples.

Signature plots
---------------

Use ``plot_signatures`` to visualize STR mutation signatures.

Each signature is shown as a bar plot of its top contributing features.

.. code-block:: python

    from str_mut_signatures import plot_signatures

    fig = plot_signatures(
        nmf_res,
        top_n=25,
    )

This is useful for:

- Interpreting biological meaning of signatures
- Comparing mutation patterns across signatures
- Selecting signatures for downstream analysis


PCA of samples
--------------

Use ``plot_pca_samples`` to see how samples cluster in PCA space based on their
exposures (or any other sample x feature matrix).

Key features:

- Takes ``NMFResult`` directly and computes PCA internally.
- Samples are colored by ``NMFResult.groups["group"]``

.. code-block:: python

    from str_mut_signatures import plot_pca_samples
    # PCA on exposures, color by "group"
    coords, var_ratio, ax = plot_pca_samples(nmf_res)

    # coords      -> DataFrame with PC1, PC2, ...
    # var_ratio   -> fraction of variance explained by each PC
    # ax          -> matplotlib Axes with the scatter plot

Why use PCA plots?

- To see whether samples separate according to known phenotypes
- To detect subgroups that may correspond to novel STR-driven biology.
- To identify outliers or potential QC issues (samples far away from the main cloud).

Using custom groups
-------------------

You can replace or augment groups manually:

.. code-block:: python

    nmf_res.groups["group"] = clinical_df.loc[nmf_res.exposures.index, "MSI_status"]

All plotting functions will automatically use the updated groups.

Command line interface
======================

Global options
--------------

- ``-v / --verbose``: Enable verbose logging.
- ``--version``: Show package version.

Extract
-------

.. code-block:: shell

    str_mut_signatures extract \
        --vcf-dir PATH \
        --out-matrix OUTPUT.tsv \
        [--ru-length] \
        [--ru {class,ru}] \
        [--ref-length] \
        [--change]

Key options:

- ``--vcf-dir``: Directory with STR-annotated, paired tumor–normal VCF files.
- ``--ru-length``: include motif length (``LEN{len(motif)}``)
- ``--ru``:

  - ``ru``: use full motif sequence.
  - ``class``: use AT-only, GC-only or mixed labeling.

- ``--ref-length``: Include reference repeat length in feature labels.
- ``--change``: Encode tumor–normal repeat-length change and restrict to somatic events.
- ``--out-matrix``: Output TSV with samples as rows and STR mutation features as columns.


Filter
------

.. code-block:: shell

    str_mut_signatures filter \
        --matrix INPUT.tsv \
        --out-matrix FILTERED.tsv \
        [--feature-method {manual,elbow,percentile}] \
        [--min-feature-total INT] \
        [--min-samples-with-feature INT] \
        [--min-sample-total INT] \
        [--feature-percentile FLOAT]

Examples:

.. code-block:: shell

    # Simple manual thresholds
    str_mut_signatures filter \
        --matrix counts_raw.tsv \
        --out-matrix counts_filtered.tsv \
        --feature-method manual \
        --min-feature-total 10 \
        --min-samples-with-feature 3 \
        --min-sample-total 0

    # Percentile-based filtering
    str_mut_signatures filter \
        --matrix counts_raw.tsv \
        --out-matrix counts_filtered.tsv \
        --feature-method percentile \
        --feature-percentile 0.9


NMF
---

.. code-block:: shell

    str_mut_signatures nmf \
        --matrix counts_filtered.tsv \
        --outdir nmf_results \
        --n-signatures 5 \
        [--max-iter 200] \
        [--random-state 0] \
        [--init nndsvd] \
        [--alpha-W 0.0] \
        [--alpha-H 0.0] \
        [--l1-ratio 0.0]

Outputs:

- ``nmf_results/signatures.tsv`` – signatures (features x K)
- ``nmf_results/exposures.tsv`` – exposures (samples x K)
- ``nmf_results/metadata.json`` – parameters and metadata


Project
-------

.. code-block:: shell

    str_mut_signatures project \
        --matrix NEW_COUNTS.tsv \
        --nmf-dir nmf_results \
        --out-exposures NEW_EXPOSURES.tsv

- ``--matrix``: New count matrix (samples x features).
- ``--nmf-dir``: Directory with an existing NMF result (``signatures.tsv``, ``metadata.json``).
- ``--out-exposures``: Output TSV with new sample exposures.


Python API
==========

Main functions
--------------

.. code-block:: python

    from str_mut_signatures import (
        parse_vcf_files,
        build_mutation_matrix,
        filter_mutation_matrix,
        run_nmf,
        save_nmf_result,
        load_nmf_result,
        project_onto_signatures,
    )

- ``parse_vcf_files(vcf_dir)`` → DataFrame of per-locus STR mutation data.
- ``build_mutation_matrix(mutations, ...)`` → samples x features count matrix.
- ``filter_mutation_matrix(matrix, ...)`` → filtered matrix + summary.
- ``run_nmf(matrix, n_signatures, ...)`` → ``NMFResult(signatures, exposures, model_params)``.
- ``save_nmf_result(result, outdir)`` / ``load_nmf_result(outdir)`` for persistence.
- ``project_onto_signatures(new_matrix, signatures, method="nnls")`` → new exposures.

Visualization helpers
---------------------

.. code-block:: python

    from str_mut_signatures import (
        plot_exposures,
        plot_pca_samples,
        plot_signatures,
    )

- ``plot_signatures(result, ...)``  
  Plot STR mutation signatures (feature loadings).

- ``plot_exposures(result, ...)``  
  Plot per-sample exposures, ordered and grouped using ``result.groups``.

- ``plot_pca_samples(result, ...)``  
  PCA of samples, colored by ``result.groups``.

Output
======

Typical outputs include:

- **Count matrices** (TSV): samples x STR mutation features.
- **Filtered matrices** (TSV): reduced feature space for robust NMF.
- **NMF signatures** and **exposures** (TSV).
- **Metadata** (JSON) describing NMF runs and parameters.
- **Figures** (from visualization helpers) for exploratory analysis and reporting.

These can be used to:

- Characterize somatic STR mutation processes.
- Compare STR signatures across cohorts.
- Associate STR signatures with clinical or genomic features.
- Apply learned STR signatures to new datasets.
- Visualize and communicate STR mutation signatures in manuscripts and presentations.


Contributing
============

Contributions are welcome!

For major changes, please open an issue first
to discuss what you’d like to change.

Please ensure:

1. All tests pass (including integration tests).
2. Code follows existing style and module structure.
3. New features include unit tests and, where appropriate, integration tests.
4. Documentation and examples are updated.


License
=======

MIT License
