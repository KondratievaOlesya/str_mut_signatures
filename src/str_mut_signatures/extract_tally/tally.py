from __future__ import annotations

from typing import Literal

import pandas as pd

RuMode = Literal[None, "length", "ru", "AT"]

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
        raise ValueError(
            "mutations_data must contain 'motif' or 'RU' column for repeat unit."
        )

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


def motif_is_at_rich(motif: str | None) -> str | pd._libs.missing.NAType:
    """
    Classify motif as AT-rich if all bases are A or T or not.

    - AT-rich: motif consists only of letters A/T.
    - non-AT-rich: motif contains any C or G (or other non-AT letters).

    Returns:
        'AT_rich' or 'non_AT_rich', or NA if motif is missing/empty.
    """
    if motif is None or pd.isna(motif):
        return pd.NA

    s = str(motif).upper()
    if not s:
        return pd.NA

    allowed = {"A", "T"}
    if all(base in allowed for base in s):
        return "AT_rich"
    return "non_AT_rich"


def make_feature(
    motif,
    ref,
    delta,
    *,
    ru: RuMode,
    ref_length: bool,
    change: bool,
):
    """
    Build a single feature key (string) for one allele / combined event.

    Parameters
    ----------
    motif
        Repeat unit (e.g. 'A', 'AT', 'AAT').
    ref
        Reference length proxy (currently from normal allele counts).
    delta
        Tumor - normal change in repeat count for this allele or combined event.
    ru
        One of:
            None    : do not include motif info
            'length': use motif length (LEN1, LEN2, ...)
            'ru'    : use full motif string (A, AT, AAT, ...)
            'AT'    : use AT-rich vs non-AT-rich classification
    ref_length
        If True, include ref in key.
    change
        If True, include non-zero delta in key and drop delta==0 events.
        If False, ignore delta for the key and do not filter by it.

    Returns
    -------
    str or NA
        Feature key like 'LEN1_10_+1', 'AT_rich_10_+2', 'A_+1', etc.,
        or NA if this event should be dropped.
    """
    if pd.isna(motif):
        return pd.NA

    parts: list[str] = []

    # Motif-based component
    if ru == "length":
        parts.append(f"LEN{len(str(motif))}")
    elif ru == "ru":
        parts.append(str(motif))
    elif ru == "AT":
        at_label = motif_is_at_rich(motif)
        if pd.isna(at_label):
            return pd.NA
        parts.append(at_label)
    elif ru is None:
        # no motif component
        pass
    else:
        raise ValueError("ru must be one of: None, 'length', 'ru', 'AT'.")

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
    ru: RuMode = "length",
    ref_length: bool = True,
    change: bool = True,
) -> pd.DataFrame:
    """
    Build somatic STR mutation count matrix from paired tumor–normal data.

    Parameters
    ----------
    mutations_data : pandas.DataFrame
        Parsed STR mutation data returned by `parse_vcf_files(...)`.

        Required columns:
            - 'sample'
            - 'normal_allele_a', 'normal_allele_b'
            - 'tumor_allele_a',  'tumor_allele_b',
            - 'motif' or 'RU'  (repeat unit sequence)
            - 'genotype_separator' (one of '|', '/', or missing)

        Phasing behaviour
        -----------------
        - If genotype_separator == '|':
            Treat genotypes as phased:
              * two allele-level events per locus (a, b).
        - Otherwise ( '/', None, missing ):
            Treat as unphased / no phasing:
              * a single **combined** event per locus
                based on total tumor vs total normal repeats.

    ru : {None, "length", "ru", "AT"}, default "length"
        Controls motif representation in feature labels:
        - None:
            Do not use motif information.
        - "length":
            Use only motif length (LEN1, LEN2, ...).
        - "ru":
            Use full motif sequence ('A', 'AT', 'AAT', ...).
        - "AT":
            Use AT-rich vs non-AT-rich classification:
                'AT_rich'    : motif consists only of A/T
                'non_AT_rich': motif contains any C/G

    ref_length : bool, default True
        If True, include a reference-length component in the feature.
        Currently this is derived from the **normal** allele counts:
          - phased: per-allele normal counts
          - unphased: combined normal count

    change : bool, default True
        If True, include tumor–normal change (delta) as part of the key
        and consider only **non-zero** changes (somatic events).
        If False, ignore delta in the key and keep all loci that pass
        basic numeric checks (presence/absence-style summaries).

    Returns
    -------
    pandas.DataFrame
        Count matrix with:
        - rows   : samples
        - columns: STR mutation feature categories
                   (defined by `ru`, `ref_length`, `change`)
        - values : counts of (allele-level or combined) events per category.
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
