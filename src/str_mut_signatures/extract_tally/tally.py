from __future__ import annotations

from typing import Literal

import pandas as pd

RuMode = Literal[None, "class", "ru"]  # none, base-class, exact RU


def validate_mutations_data(df: pd.DataFrame) -> tuple[str, bool]:
    """
    Validate the input DataFrame and return:
      - motif column name ('motif' or 'RU')
      - whether genotype_separator column is present
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("mutations_data must be a pandas.DataFrame")

    required_cols = {
        "sample",
        "normal_allele_a",
        "normal_allele_b",
        "tumor_allele_a",
        "tumor_allele_b",
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"mutations_data is missing required columns: {missing}")

    if "motif" in df.columns:
        motif_col = "motif"
    elif "RU" in df.columns:
        motif_col = "RU"
    else:
        raise ValueError("mutations_data must contain 'motif' or 'RU' column for repeat unit.")

    has_genotype_sep = "genotype_separator" in df.columns

    return motif_col, has_genotype_sep


def is_phased(genotype_separator: str | None) -> bool:
    """
    Return True if genotype separator indicates phased genotypes ('|').
    """
    return genotype_separator == "|"


def compute_changes_for_row(row: pd.Series) -> pd.Series:
    """
    Compute allele-level or combined tumor–normal changes for a single row.

    - If phased (GT uses '|'):
        change_a = tumor_allele_a - normal_allele_a
        change_b = tumor_allele_b - normal_allele_b
        ref_a    = normal_allele_a
        ref_b    = normal_allele_b

    - If unphased or no phasing info:
        We only track a single combined change:
            total_normal = normal_allele_a + normal_allele_b
            total_tumor  = tumor_allele_a + tumor_allele_b
            change_total = total_tumor - total_normal

        This is stored in:
            change_a = change_total
            ref_a    = total_normal
            change_b = NA
            ref_b    = NA
    """
    # Extract genotype separator if present
    genotype_separator = row.get("genotype_separator", None)
    phased = is_phased(genotype_separator)

    try:
        n_a = int(row["normal_allele_a"])
        n_b = int(row["normal_allele_b"])
        t_a = int(row["tumor_allele_a"])
        t_b = int(row["tumor_allele_b"])
    except Exception:
        # If parsing fails, treat as missing
        return pd.Series(
            {
                "change_a": pd.NA,
                "change_b": pd.NA,
                "ref_a": pd.NA,
                "ref_b": pd.NA,
            }
        )

    if phased:
        # Allele-specific changes
        change_a = t_a - n_a
        change_b = t_b - n_b
        ref_a = n_a
        ref_b = n_b
    else:
        # Combined change only
        total_normal = n_a + n_b
        total_tumor = t_a + t_b
        change_total = total_tumor - total_normal

        change_a = change_total
        ref_a = total_normal
        change_b = pd.NA
        ref_b = pd.NA

    return pd.Series(
        {
            "change_a": change_a,
            "change_b": change_b,
            "ref_a": ref_a,
            "ref_b": ref_b,
        }
    )


def motif_base_class(motif: str | None) -> str | pd._libs.missing.NAType:
    """
    Classify a repeat-unit motif by nucleotide base composition.

    The motif is assigned to one of three base-composition classes based on
    the presence of A/T and G/C nucleotides.

    Parameters
    ----------
    motif : str or None
        Repeat-unit sequence (e.g. ``"A"``, ``"AT"``, ``"AAT"``).

    Returns
    -------
    str or pandas.NA
        Base-composition class of the motif:

        - ``"AT_only"`` : motif contains only A/T bases
        - ``"GC_only"`` : motif contains only G/C bases
        - ``"mixed"`` : motif contains both A/T and G/C bases
        Returns ``pandas.NA`` if the motif is missing, empty, or contains
        non-ACGT characters.
    """
    if motif is None or pd.isna(motif):
        return pd.NA

    s = str(motif).strip().upper()
    if not s:
        return pd.NA

    allowed = set("ACGT")
    chars = set(s)

    # handle invalid characters
    if not chars.issubset(allowed):
        return pd.NA

    at = {"A", "T"}
    gc = {"G", "C"}

    has_at = len(chars & at) > 0
    has_gc = len(chars & gc) > 0

    if has_at and not has_gc:
        return "AT_only"
    if has_gc and not has_at:
        return "GC_only"
    return "mixed"


def make_feature(
    motif,
    ref,
    delta,
    *,
    ru_length: bool,
    ru: RuMode,
    ref_length: bool,
    change: bool,
):
    """
    Construct a single STR mutation feature key for one allele or combined event.

    The feature key is composed of optional components describing repeat-unit
    length, repeat-unit content, reference length, and somatic change. Components
    are concatenated using underscores.

    Parameters
    ----------
    motif : str or pandas.NA
        Repeat unit sequence (e.g. ``"A"``, ``"AT"``, ``"AAT"``).

    ref : int or pandas.NA
        Reference repeat length, typically derived from the normal allele
        repeat count.

    delta : int or pandas.NA
        Tumor–normal change in repeat count for this allele or combined event.

    ru_length : bool
        If True, include the repeat-unit length as ``LEN{len(motif)}``
        in the feature key.

    ru : {None, "ru", "class"}
        Controls how repeat-unit *content* is represented in the feature key.

        - ``None`` :
          Do not include repeat-unit content.
        - ``"ru"`` :
          Include the full repeat-unit sequence (e.g. ``"A"``, ``"AT"``).
        - ``"class"`` :
          Include the base-composition class of the repeat unit:

          - ``AT_only`` : motif contains only A/T
          - ``GC_only`` : motif contains only G/C
          - ``mixed`` : mixed A/T and G/C

    ref_length : bool
        If True, include the reference repeat length in the feature key.

    change : bool
        If True, include the tumor–normal repeat count change (delta) in the
        feature key and discard events with ``delta == 0``.

        If False, ignore delta and retain all loci that pass basic numeric
        checks.

    Returns
    -------
    str or pandas.NA
        Feature key string (e.g. ``"LEN1_AT_only_10_+1"``, ``"GC_only_12_-2"``,
        ``"A_+1"``), or ``pandas.NA`` if the event should be discarded.
    """
    if pd.isna(motif):
        return pd.NA

    parts: list[str] = []
    motif_s = str(motif).strip().upper()
    if not motif_s:
        return pd.NA

    # RU length component (optional)
    if ru_length:
        parts.append(f"LEN{len(motif_s)}")

    # RU content component (optional)
    if ru == "class":
        ru_cls = motif_base_class(motif_s)  # AT_only / GC_only / mixed
        if pd.isna(ru_cls):
            return pd.NA
        parts.append(ru_cls)
    elif ru == "ru":
        parts.append(motif_s)
    elif ru is not None:
        raise ValueError("ru must be one of: None, 'class', 'ru'.")

    # Reference length component
    if ref_length:
        if pd.isna(ref):
            return pd.NA
        parts.append(str(int(ref)))

    # Somatic change component
    if change:
        if pd.isna(delta):
            return pd.NA
        d = int(delta)
        # Only count true somatic events
        if d == 0:
            return pd.NA
        sign = "+" if d > 0 else ""
        parts.append(f"{sign}{d}")

    # If everything was turned off (no ru, no ref_length, no change)
    if not parts:
        return pd.NA

    return "_".join(parts)


def build_mutation_matrix(
    mutations_data: pd.DataFrame,
    *,
    ru_length: bool = True,
    ru: RuMode = None,
    ref_length: bool = True,
    change: bool = True,
) -> pd.DataFrame:
    """
    Build a somatic STR mutation count matrix from paired tumor–normal data.

    This function converts per-locus STR mutation calls into a sample-by-feature
    count matrix. Feature definitions are controlled by repeat-unit length,
    repeat-unit content, reference length, and somatic change options.

    Parameters
    ----------
    mutations_data : pandas.DataFrame
        Parsed STR mutation data, typically returned by
        :func:`parse_vcf_files`.

        Required columns include:

        - ``sample``
        - ``normal_allele_a``, ``normal_allele_b``
        - ``tumor_allele_a``, ``tumor_allele_b``
        - ``motif`` or ``RU`` (repeat unit sequence)
        - ``genotype_separator`` (``'|'``, ``'/'``, or missing)

    ru_length : bool, default True
        If True, include the repeat-unit length as ``LEN{len(motif)}``
        in the feature key.

    ru : {None, "class", "ru"}, default None
        Controls how repeat-unit *content* is represented in the feature key.

        - ``None`` :
          Do not include repeat-unit content.
        - ``"ru"`` :
          Use the full repeat-unit sequence (e.g. ``A``, ``AT``, ``AAT``).
        - ``"class"`` :
          Use base-composition class of the repeat unit:

          - ``AT_only`` : motif contains only A/T
          - ``GC_only`` : motif contains only G/C
          - ``mixed`` : mixed A/T and G/C

    ref_length : bool, default True
        If True, include a reference-length component derived from the
        normal allele repeat counts.

        - Phased genotypes: per-allele normal repeat count
        - Unphased genotypes: combined normal repeat count

    change : bool, default True
        If True, include the tumor–normal repeat count change (delta) in
        the feature key and retain only non-zero changes (somatic events).

        If False, ignore delta and retain all loci that pass basic numeric
        checks, producing presence/absence-style summaries.

    Returns
    -------
    pandas.DataFrame
        STR mutation count matrix with:

        - rows: samples
        - columns: STR mutation feature categories
        - values: counts of allele-level or combined STR mutation events

    Notes
    -----
    Phasing behavior is determined by ``genotype_separator``:

    - ``'|'`` :
      Genotypes are treated as phased, producing two allele-level events
      per locus.
    - ``'/'`` or missing :
      Genotypes are treated as unphased, producing a single combined event
      per locus based on total tumor vs. normal repeat counts.
    """
    df = mutations_data.copy()
    motif_col, has_genotype_sep = validate_mutations_data(df)

    # Compute allele-level / combined changes and reference lengths
    changes = df.apply(compute_changes_for_row, axis=1)
    df[["change_a", "change_b", "ref_a", "ref_b"]] = changes

    # Build feature labels for each allele / combined event
    df["mutation_type_a"] = df.apply(
        lambda row: make_feature(
            motif=row[motif_col],
            ref=row["ref_a"],
            delta=row["change_a"],
            ru_length=ru_length,
            ru=ru,
            ref_length=ref_length,
            change=change,
        ),
        axis=1,
    )

    df["mutation_type_b"] = df.apply(
        lambda row: make_feature(
            motif=row[motif_col],
            ref=row["ref_b"],
            delta=row["change_b"],
            ru_length=ru_length,
            ru=ru,
            ref_length=ref_length,
            change=change,
        ),
        axis=1,
    )

    # Long format: one row per (sample, allele/combo-level mutation_type)
    df_long = pd.melt(
        df,
        id_vars=["sample"],
        value_vars=["mutation_type_a", "mutation_type_b"],
        var_name="allele_type",
        value_name="mutation_type",
    )

    # Drop entries without a valid feature (e.g. non-somatic when change=True)
    df_long = df_long.dropna(subset=["mutation_type"])

    # If nothing left (e.g. no somatic events), return empty matrix
    if df_long.empty:
        return pd.DataFrame()

    # Count matrix: samples x mutation_type
    mutation_counts = (
        df_long.groupby(["sample", "mutation_type"])
        .size()
        .unstack(fill_value=0)
        .sort_index(axis=0)
        .sort_index(axis=1)
    )

    return mutation_counts
