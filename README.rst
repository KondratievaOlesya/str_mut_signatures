==================
str_mut_signatures
==================

STR Mutation Signature Analysis.

Python package for analysis of Short Tandem Repeat (STR) mutation signatures from VCF files.
It extracts somatic STR mutation events from paired tumor–normal VCFs, builds
count matrices, and performs NMF-based signature decomposition and visualization.

Contents
========

- Installation
- Quick start
- Input format
- Somatic STR calls (tumor–normal)
- Annotating standard VCFs
- Matrix construction
- Command line interface
- Python API
- Output
- Contributions
- License


Installation
============

From PyPI
--------

Package is available through `PyPI <https://pypi.org/project/str_mut_signatures/>`_. To install, run:

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

Command Line
------------

1. Extract somatic STR mutations from paired tumor–normal VCFs and build a count matrix:

.. code-block:: shell

    str_mut_signatures extract \
        --vcf-dir data/vcfs/ \
        --out-matrix counts_len1.tsv \
        --ru length \
        --ref-length \
        --change

This produces features of the form:

.. code-block:: text

    LEN{motif_length}_{ref_length}_{change}

e.g. ``LEN1_10_+1`` for a 1-bp motif, reference length 10, and +1 repeat unit change
in tumor vs normal.

2. Run NMF decomposition:

.. code-block:: shell

    str_mut_signatures nmf \
        --matrix counts_len1.tsv \
        --outdir nmf_results


Python Library Usage
--------------------

.. code-block:: python

    from str_mut_signatures import (
        parse_vcf_files,
        build_mutation_matrix,
        run_nmf,
    )

    # Parse annotated paired tumor–normal VCF files
    mutations = parse_vcf_files("vcf_directory/")

    # Build a mutation count matrix:
    # ru:
    #   None   -> ignore motif
    #   "length" -> use only motif length (LEN1, LEN2, ...)
    #   "ru"   -> use full repeat unit sequence (e.g. AT, AAT)
    # ref_length:
    #   include reference repeat length as a feature component
    # change:
    #   include tumor–normal repeat-length change
    matrix = build_mutation_matrix(
        mutations,
        ru="length",
        ref_length=True,
        change=True,
    )

    # Example resulting column name:
    # LEN1_{ref_length}_{change}
    # e.g. LEN1_10_+1

    matrix.to_csv("counts_matrix.tsv", sep="\t")

    # Run NMF
    W, H = run_nmf(matrix, n_signatures=5)


Input format
============

To be processed by ``str_mut_signatures``, VCF files must:

1. Contain **paired samples** (normal and tumor) per record.
2. Be annotated with STR-specific fields that describe the repeat unit and allele-level
   repeat counts.

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
  - ``REF``: Reference repeat count (copy number of the motif in the reference genome)
    or equivalent information if available.

  **FORMAT fields**

  - ``REPCN``: Comma-separated repeat copy numbers for each allele in that sample.

Example schema
--------------

.. code-block:: text

    ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
    ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference repeat count">
    ##FORMAT=<ID=REPCN,Number=R,Type=Integer,Description="Per-allele repeat copy number">
    #CHROM POS  ID REF ALT QUAL FILTER INFO        FORMAT   NORMAL   TUMOR
    chr1   100 .  A   AT  .    .     RU=A;REF=10  REPCN    10,10    10,11

From this, ``str_mut_signatures``:

- Compares NORMAL vs TUMOR ``REPCN``.
- Identifies loci where tumor repeat copy number differs from normal.
- Encodes the net repeat-length **change** as tumor–normal (e.g. ``+1``).
- Uses only these somatic events for downstream count matrices and signatures.


Somatic STR calls (tumor–normal)
================================

Key points:

- The package focuses on **somatic** STR mutations.
- For each locus, tumor and normal alleles are compared:
  
  - If there is no difference between tumor and normal (based on ``REPCN``), the
    site is ignored.
  - If there is a difference, a somatic STR event is recorded.

- The ``change`` feature encodes **tumor–normal** repeat-length difference, not
  reference–sample:

  .. code-block:: text

      change = f(REPCN_tumor, REPCN_normal)

  Implementation details (e.g. handling of heterozygous states) follow the
  library’s internal definition; the core idea is that only tumor–normal
  differences contribute to the matrix.


Annotating standard VCFs
========================

If your VCFs lack ``RU``, ``REF``, or ``REPCN``, you can annotate them using
the companion tool ``strvcf_annotator``:

- Takes standard VCF + STR reference.
- Produces STR-annotated VCFs compatible with ``str_mut_signatures``.

For details see: ``strvcf_annotator``_.

.. _strvcf_annotator: https://github.com/acg-team/strvcf_annotator


Matrix construction
===================

``build_mutation_matrix`` provides a flexible way to define feature space
(columns) using simple flags.

Core components
---------------

Given:

- ``RU``: repeat unit sequence
- ``len(RU)``: motif length
- ``REF``: reference repeat count
- ``change``: tumor–normal repeat-length change at that locus

you can select:

- ``ru``:
  
  - ``None``:
    
    - Do not use motif information.

  - ``"length"``:
    
    - Use only motif length.
    - Features start with ``LEN{motif_length}``.
    - Example: ``LEN1_10_+1`` for motif length 1, REF=10, change=+1.

  - ``"ru"``:
    
    - Use full repeat unit sequence.
    - Example: ``A_10_+1``, ``AT_20_-2``.

- ``ref_length`` (bool):

  - If ``True``, include the reference repeat length as part of the feature key.

- ``change`` (bool):

  - If ``True``, include tumor–normal change as part of the feature key.
  - Only somatic events (non-zero change) are counted.

Examples
--------

1. Motif length + ref length + somatic change:

.. code-block:: python

    m = build_mutation_matrix(
        mutations,
        ru="length",
        ref_length=True,
        change=True,
    )
    # Columns: LEN{motif_length}_{ref_length}_{change}
    # e.g. LEN1_10_+1

2. Full motif + change only:

.. code-block:: python

    m = build_mutation_matrix(
        mutations,
        ru="ru",
        ref_length=False,
        change=True,
    )
    # Columns: {RU}_{change}
    # e.g. AT_+2

3. Motif length only (no change, e.g. for presence/absence-style summaries):

.. code-block:: python

    m = build_mutation_matrix(
        mutations,
        ru="length",
        ref_length=False,
        change=False,
    )
    # Columns: LEN1, LEN2, ...


Command line interface
======================

Extract
-------

.. code-block:: shell

    str_mut_signatures extract \
        --vcf-dir PATH \
        --out-matrix OUTPUT.tsv \
        [--ru {none,length,ru}] \
        [--ref-length] \
        [--change]

Key options:

- ``--vcf-dir``: Directory with STR-annotated, paired tumor–normal VCF files.
- ``--ru``:

  - ``none``: ignore motif.
  - ``length``: use motif length (LEN1, LEN2, ...).
  - ``ru``: use full motif sequence.

- ``--ref-length``: Include reference repeat length in feature labels.
- ``--change``: Encode tumor–normal repeat-length change and restrict to somatic events.
- ``--out-matrix``: Output TSV with samples as rows and STR mutation features as columns.


NMF
---

.. code-block:: shell

    str_mut_signatures nmf \
        --matrix counts.tsv \
        --outdir nmf_results \
        --n-signatures 5

Outputs:

- Signature profiles (W)
- Sample exposures (H)
- Basic diagnostic plots


Python API
==========

.. code-block:: python

    from str_mut_signatures import parse_vcf_files, build_mutation_matrix, run_nmf

    mutations = parse_vcf_files("vcf_directory/")

    m = build_mutation_matrix(
        mutations,
        ru="length",
        ref_length=True,
        change=True,
    )

    W, H = run_nmf(m, n_signatures=5)


Output
======

- Count matrices (TSV): samples x STR mutation features.
- NMF signatures and exposures.
- Visualization plots (signatures, exposures).
- Basic analysis metrics.

These can be used to:

- Characterize somatic STR mutation processes.
- Compare STR signatures across cohorts.
- Associate STR signatures with clinical or genomic features.


Documentation
==============

* `API Documentation <docs/API.md>`_
* `Examples <examples/>`_

Contributing
============

Contributions are welcome! 
For major changes, please open an issue first 
to discuss what you’d like to change.
Please ensure:

1. All tests pass
2. Code follows existing style
3. New features include tests
4. Documentation is updated

License
============

MIT License