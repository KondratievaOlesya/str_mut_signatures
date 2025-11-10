import pandas as pd


def build_mutation_matrix(
    mutations_data,
    ru="length",
    ref_length=True,
    change=True,
):
    """
    Build somatic STR mutation count matrix from paired tumor–normal VCF data.

    Parameters
    ----------
    mutations_data : dict
        Parsed STR mutation data returned by `parse_vcf_files(...)`.
        Expected to originate from STR-annotated VCFs with:
        - paired samples per record, where:
            - the FIRST sample is NORMAL
            - the SECOND sample is TUMOR
        - INFO fields:
            - RU   : repeat unit sequence (e.g. "A", "AT", "AAT")
            - REF  : reference repeat count
        - FORMAT fields:
            - REPCN: per-allele repeat copy numbers for each sample
        Only loci with differences between tumor and normal (based on REPCN-derived
        repeat lengths) are considered; i.e. the matrix is built from SOMATIC STR mutations.

    ru : {None, "length", "ru"}, default "length"
        Controls how the repeat unit is represented in feature labels:
        - None:
            Do not use motif information.
        - "length":
            Use only the motif length.
            Features start with `LEN{motif_length}`.
        - "ru":
            Use full repeat unit sequence.
            Features start with the motif itself (e.g. `A`, `AT`, `AAT`).

    ref_length : bool, default True
        If True, include the reference repeat length (REF-derived length) as part
        of the feature key.

    change : bool, default True
        If True, include the TUMOR–NORMAL repeat-length change as part of the feature key.
        The change is computed from tumor vs normal allele copy numbers, **not**
        from reference vs sample. Only non-zero tumor–normal differences are counted.

        Example (with ru="length", ref_length=True, change=True):
            Feature label format:
                LEN{motif_length}_{ref_length}_{change}
            e.g. `LEN1_10_+1` for:
                - motif length = 1
                - reference length = 10
                - tumor has +1 repeat unit relative to normal.

    Returns
    -------
    pandas.DataFrame
        Count matrix with:
        - rows: samples
        - columns: STR mutation feature categories defined by (ru, ref_length, change)
        - entries: counts of somatic STR mutations per category.
    """
    # --- Basic checks & setup -------------------------------------------------
    if not isinstance(mutations_data, pd.DataFrame):
        raise TypeError("mutations_data must be a pandas.DataFrame")

    required_cols = {
        "sample",
        "normal_allele_a",
        "normal_allele_b",
        "tumor_allele_a",
        "tumor_allele_b",
    }
    missing = required_cols - set(mutations_data.columns)
    if missing:
        raise ValueError(f"mutations_data is missing required columns: {missing}")

    # Motif column: allow either 'motif' or 'RU'
    if "motif" in mutations_data.columns:
        motif_col = "motif"
    elif "RU" in mutations_data.columns:
        motif_col = "RU"
    else:
        raise ValueError("mutations_data must contain 'motif' or 'RU' column for repeat unit.")

    if ru not in (None, "length", "ru"):
        raise ValueError("ru must be one of: None, 'length', 'ru'.")

    df = mutations_data.copy()

    # --- Helper: compute per-allele tumor–normal changes ----------------------
    def compute_allele_changes(row):
        try:
            tumor = sorted([
                int(row["tumor_allele_a"]),
                int(row["tumor_allele_b"]),
            ])
            normal = sorted([
                int(row["normal_allele_a"]),
                int(row["normal_allele_b"]),
            ])
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

        change_a = tumor[0] - normal[0]
        change_b = tumor[1] - normal[1]
        ref_a = normal[0]
        ref_b = normal[1]

        return pd.Series(
            {
                "change_a": change_a,
                "change_b": change_b,
                "ref_a": ref_a,
                "ref_b": ref_b,
            }
        )

    df[["change_a", "change_b", "ref_a", "ref_b"]] = df.apply(
        compute_allele_changes,
        axis=1,
    )

    # --- Helper: build feature label for a single allele ----------------------
    def make_feature(motif, ref, delta):
        # If any critical piece is missing, skip
        if pd.isna(motif):
            return pd.NA

        # Handle motif representation
        parts = []

        if ru == "length":
            parts.append(f"LEN{len(str(motif))}")
        elif ru == "ru":
            parts.append(str(motif))
        elif ru is None:
            # no motif info
            pass

        # Reference length (from normal)
        if ref_length:
            if pd.isna(ref):
                return pd.NA
            parts.append(str(int(ref)))

        # Somatic change (tumor - normal)
        if change:
            # Only keep true somatic events
            if pd.isna(delta) or int(delta) == 0:
                return pd.NA
            d = int(delta)
            sign = "+" if d > 0 else ""
            parts.append(f"{sign}{d}")

        # If user turned everything off (no ru, no ref_length, no change)
        if not parts:
            return pd.NA

        return "_".join(parts)

    # --- Construct mutation types per allele ---------------------------------
    df["mutation_type_a"] = df.apply(
        lambda row: make_feature(
            row[motif_col],
            row["ref_a"],
            row["change_a"],
        ),
        axis=1,
    )

    df["mutation_type_b"] = df.apply(
        lambda row: make_feature(
            row[motif_col],
            row["ref_b"],
            row["change_b"],
        ),
        axis=1,
    )

    # --- Long format: one row per (sample, allele-level mutation_type) -------
    df_long = pd.melt(
        df,
        id_vars=["sample"],
        value_vars=["mutation_type_a", "mutation_type_b"],
        var_name="allele_type",
        value_name="mutation_type",
    )

    # Drop entries with no feature (e.g. non-somatic when change=True)
    df_long = df_long.dropna(subset=["mutation_type"])

    # If nothing left (e.g. no somatic events), return empty matrix
    if df_long.empty:
        return pd.DataFrame()

    # --- Count matrix: samples x mutation_type --------------------------------
    mutation_counts = (
        df_long
        .groupby(["sample", "mutation_type"])
        .size()
        .unstack(fill_value=0)
        .sort_index(axis=0)
        .sort_index(axis=1)
    )

    return mutation_counts
