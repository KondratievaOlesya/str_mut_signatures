from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd

RuMode = Literal[None, "class", "ru"]  # none, base-class, exact RU


def validate_mutations_data(df: pd.DataFrame) -> str:
    """
    Validate the input DataFrame and return:
      - motif column name ('motif' or 'RU')
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

    return motif_col



def motif_base_class(motif: str | None) -> str | None:
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
        - Unphased genotypes: sorts allele pairs and compares them element-wise

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
      Genotypes are treated as unphased, sorts allele pairs and compares them element-wise
    """
    df = mutations_data.copy()
    motif_col = validate_mutations_data(df)
    # Convert allele columns to numeric once
    allele_cols = [
        "normal_allele_a",
        "normal_allele_b",
        "tumor_allele_a",
        "tumor_allele_b",
    ]
    for col in allele_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=allele_cols).copy()
    if df.empty:
        return pd.DataFrame()

    for col in allele_cols:
        df[col] = df[col].astype(int)

    # Create phased mask once
    if "genotype_separator" in df.columns:
        phased_mask = df["genotype_separator"].eq("|")
    else:
        phased_mask = pd.Series(False, index=df.index)

    # Prepare output columns
    df["change_a"] = pd.NA
    df["change_b"] = pd.NA
    df["ref_a"] = pd.NA
    df["ref_b"] = pd.NA

    # --- Phased rows: compare alleles by position ---
    phased_idx = df.index[phased_mask]
    if len(phased_idx) > 0:
        df.loc[phased_idx, "change_a"] = (
            df.loc[phased_idx, "tumor_allele_a"] - df.loc[phased_idx, "normal_allele_a"]
        )
        df.loc[phased_idx, "change_b"] = (
            df.loc[phased_idx, "tumor_allele_b"] - df.loc[phased_idx, "normal_allele_b"]
        )
        df.loc[phased_idx, "ref_a"] = df.loc[phased_idx, "normal_allele_a"]
        df.loc[phased_idx, "ref_b"] = df.loc[phased_idx, "normal_allele_b"]

    # --- Unphased rows: sort allele pairs before comparison ---
    unphased_idx = df.index[~phased_mask]
    if len(unphased_idx) > 0:
        normal_pairs = np.sort(
            df.loc[unphased_idx, ["normal_allele_a", "normal_allele_b"]].to_numpy(),
            axis=1,
        )
        tumor_pairs = np.sort(
            df.loc[unphased_idx, ["tumor_allele_a", "tumor_allele_b"]].to_numpy(),
            axis=1,
        )

        df.loc[unphased_idx, "ref_a"] = normal_pairs[:, 0]
        df.loc[unphased_idx, "ref_b"] = normal_pairs[:, 1]
        df.loc[unphased_idx, "change_a"] = tumor_pairs[:, 0] - normal_pairs[:, 0]
        df.loc[unphased_idx, "change_b"] = tumor_pairs[:, 1] - normal_pairs[:, 1]
    # Use nullable integer dtype
    for col in ["change_a", "change_b", "ref_a", "ref_b"]:
        df[col] = pd.array(df[col], dtype="Int64")

    df = df.loc[~(df["change_a"].eq(0) & df["change_b"].eq(0))].copy()
    # Fallback: no features selected -> just count mutations per sample
    if not ru_length and ru is None and not ref_length and not change:
        df["event_a"] = df["change_a"].fillna(0).ne(0).astype(int)
        df["event_b"] = df["change_b"].fillna(0).ne(0).astype(int)

        return (
            df.assign(mutation_count=df["event_a"] + df["event_b"])
            .groupby("sample", sort=True)["mutation_count"]
            .sum()
            .to_frame()
        )

    # Build motif-derived components once
    motif = df[motif_col].astype("string").str.strip().str.upper()

    valid_motif = motif.notna() & motif.ne("")
    if ru == "class":
        motif_component = motif.map(motif_base_class)
        valid_motif &= motif_component.notna()
    elif ru == "ru":
        motif_component = motif
        valid_motif &= motif.notna()
    elif ru is None:
        motif_component = pd.Series(pd.NA, index=df.index, dtype="object")
    else:
        raise ValueError("ru must be one of: None, 'class', 'ru'.")

    if ru_length:
        motif_len_component = "LEN" + motif.str.len().astype("Int64").astype("string")
        valid_motif &= motif.notna()
    else:
        motif_len_component = pd.Series(pd.NA, index=df.index, dtype="object")

    # Create long format first, then build features vectorized
    left = df[["sample"]].copy()
    left["motif_valid"] = valid_motif
    left["motif_component"] = motif_component
    left["motif_len_component"] = motif_len_component
    left["ref"] = df["ref_a"]
    left["delta"] = df["change_a"]

    right = df[["sample"]].copy()
    right["motif_valid"] = valid_motif
    right["motif_component"] = motif_component
    right["motif_len_component"] = motif_len_component
    right["ref"] = df["ref_b"]
    right["delta"] = df["change_b"]

    df_long = pd.concat([left, right], axis=0, ignore_index=True)

    # Validity mask
    keep = df_long["motif_valid"].fillna(False)

    if ref_length:
        keep &= df_long["ref"].notna()

    if change:
        keep &= df_long["delta"].notna()
        keep &= df_long["delta"].ne(0)

    df_long = df_long.loc[keep].copy()

    # If no valid features remain, return empty DataFrame
    if df_long.empty:
        return pd.DataFrame()
    # Build feature strings from enabled components
    parts = []

    if ru_length:
        parts.append(df_long["motif_len_component"].astype("string"))

    if ru in {"class", "ru"}:
        parts.append(df_long["motif_component"].astype("string"))

    if ref_length:
        parts.append(df_long["ref"].astype("Int64").astype("string"))

    if change:
        delta_int = df_long["delta"].astype("Int64")
        delta_str = delta_int.astype(str)
        delta_str = np.where(delta_int > 0, "+" + delta_str, delta_str)
        parts.append(pd.Series(delta_str, index=df_long.index, dtype="string"))

    if not parts:
        raise ValueError(
            "At least one feature component must be enabled: "
            "ru_length, ru, ref_length, or change."
        )

    feature = parts[0]
    for p in parts[1:]:
        feature = feature + "_" + p

    df_long["mutation_type"] = feature

    mutation_counts = (
        df_long.groupby(["sample", "mutation_type"], sort=True)
        .size()
        .unstack(fill_value=0)
        .sort_index(axis=0)
        .sort_index(axis=1)
    )

    return mutation_counts
