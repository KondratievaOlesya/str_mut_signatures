from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
from trtools.utils.utils import GetCanonicalMotif

from .validate import validate_vcf


def normalize_motif(raw_motif: str | None) -> str:
    """
    Normalize motif string.

    Normalization steps:

    - Convert to uppercase
    - Remove spaces and non-letter characters
    - Canonicalize using ``trtools.utils.utils.GetCanonicalMotif``

    Returns an empty string if nothing usable remains.

    Parameters
    ----------
    raw_motif : str or None
        Input motif string to normalize. If ``None``, an empty string is returned.

    Returns
    -------
    str
        Normalized canonical motif, or an empty string if invalid.
    """
    if raw_motif is None or pd.isna(raw_motif):
        return ""

    # Uppercase and keep only letters (A–Z)
    cleaned = "".join(ch for ch in str(raw_motif).upper() if ch.isalpha())

    if not cleaned:
        return ""

    # Canonical motif
    canonical = GetCanonicalMotif(cleaned)
    return canonical


def parse_info(info_field: str) -> dict:
    info = {}
    for item in info_field.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            info[key] = val
    return info


def parse_copy_number(cn_str: str):
    """
    Parse copy-number field (``REPCN`` / ``REPLEN``) into two alleles.

    The input is expected to be a comma-separated string, e.g. ``"10,11"``.

    Rules:

    - If a single value is provided, it is treated as homozygous (``a == b``).
    - If two values are provided, they are returned as-is.
    - If more than two values are provided, returns ``(".", ".")`` (caller can skip).

    Parameters
    ----------
    cn_str : str
        Copy-number string to parse (e.g. ``"10,11"``).

    Returns
    -------
    tuple[str, str]
        Two allele values as strings.
    """
    parts = cn_str.split(",")
    if len(parts) == 2:
        return parts[0], parts[1]
    elif len(parts) == 1:
        return parts[0], parts[0]
    else:
        return ".", "."


def open_maybe_gzip(path: str | Path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def process_vcf_to_rows(
    path: str | Path,
    *,
    filter_by_pass: bool = True,
    filter_by_perfect: bool = True,
):
    """
    Parse a single STR-annotated VCF into row dictionaries.

    Supports:

    - GangSTR: uses ``FORMAT/REPCN`` as copy number
    - conSTRain: uses ``FORMAT/REPLEN`` as copy number
    - VCF annotated with ``strvcf_annotator`` (``INFO/RU``, ``INFO/REF``, ``FORMAT/REPCN``)

    Filtering options
    -----------------
    filter_by_pass
        If ``True`` (default), keep only records with ``FILTER == "PASS"``.
        If ``False``, ignore the ``FILTER`` field.
    filter_by_perfect
        If ``True`` (default), and ``INFO/PERFECT`` is present, keep only records
        where ``PERFECT != "FALSE"`` (i.e. skip variants where ``PERFECT == "FALSE"``).
        If ``False``, ignore the ``PERFECT`` flag completely.

    Assumptions after validation
    ----------------------------
    - First sample column after ``FORMAT`` is NORMAL (index 9 in standard VCF).
    - Second sample column after ``FORMAT`` is TUMOR (index 10).
    - STR annotations are present in ``INFO`` / ``FORMAT``.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the STR-annotated VCF file.
    filter_by_pass : bool, optional
        Whether to keep only records with ``FILTER == "PASS"``.
    filter_by_perfect : bool, optional
        Whether to filter by ``INFO/PERFECT`` when present.

    Returns
    -------
    list[dict]
        List of dictionaries, one per parsed STR record.
    """
    path = Path(path)

    # --- Validate VCF first ---
    vres = validate_vcf(path)

    if not vres.has_str_annotations:
        raise ValueError(
            f"VCF '{path}' is missing required STR annotations "
            f"(INFO/RU, INFO/REF, and either FORMAT/REPCN or FORMAT/REPLEN). "
            f"Annotate it with 'strvcf_annotator' or use a supported STR caller "
            f"(GangSTR / conSTRain) before using str_mut_signatures."
        )

    if vres.copy_number_field is None:
        # Should not happen if has_str_annotations is True, but keep it explicit
        raise ValueError(
            f"VCF '{path}' does not define a supported copy-number FORMAT field "
            f"(expected REPCN for GangSTR or REPLEN for conSTRain)."
        )

    cn_field = vres.copy_number_field  # "REPCN" or "REPLEN"

    if not vres.has_paired_samples:
        raise ValueError(
            f"VCF '{path}' must contain at least two samples "
            f"(normal and tumor) after the FORMAT column."
        )

    if vres.genotype_separator is None:
        raise ValueError(
            f"VCF '{path}' does not have a usable GT field in the first record "
            f"(could not determine phased '|' vs unphased '/')."
        )

    genotype_separator = vres.genotype_separator  # '|' or '/'

    rows = []
    filter_not_passed_count = 0
    non_perfect_count = 0
    written_variants = 0

    # We are using filename as sample identifier
    sample_name = path.name.replace(".vcf", "").replace(".vcf.gz", "")

    header_cols = None
    normal_idx = None
    tumor_idx = None
    format_fields = None

    with open_maybe_gzip(path) as f:
        for line in f:
            # Header with column names
            if line.startswith("#CHROM"):
                header_cols = line.strip().split()
                if len(header_cols) < 11:
                    raise ValueError(
                        f"VCF '{path}' must have at least 2 samples (normal, tumor) "
                        f"in the header after FORMAT."
                    )
                # By convention: first = NORMAL, second = TUMOR
                normal_idx = 9
                tumor_idx = 10
                continue

            # Other header lines
            if line.startswith("#"):
                continue

            # Data line
            cols = line.strip().split()
            if len(cols) < 10:
                # malformed, skip
                continue

            # FILTER == PASS (optional)
            if filter_by_pass and cols[6] != "PASS":
                filter_not_passed_count += 1
                continue

            info = parse_info(cols[7])

            # PERFECT calls only (optional, only applies if PERFECT exists)
            if filter_by_perfect and info.get("PERFECT", "") == "FALSE":
                non_perfect_count += 1
                continue

            chrom = cols[0]
            pos = cols[1]
            tmp_id = f"{chrom}_{pos}"

            end = info.get("END", "")
            period = info.get("PERIOD", "")
            ref = info.get("REF", "")
            motif = info.get("RU", "")
            motif = normalize_motif(motif)
            # FORMAT & samples
            if format_fields is None:
                format_fields = cols[8].split(":")

            try:
                normal_sample = cols[normal_idx]
                tumor_sample = cols[tumor_idx]
            except IndexError as err:
                raise ValueError(
                    f"VCF '{path}' does not contain expected normal/tumor sample "
                    f"columns at indices 9 and 10."
                ) from err

            # parse NORMAL copy numbers
            normal_values = normal_sample.split(":")
            normal_fmt = dict(zip(format_fields, normal_values))
            n_a, n_b = parse_copy_number(normal_fmt.get(cn_field, ".,."))

            # parse TUMOR copy numbers
            tumor_values = tumor_sample.split(":")
            tumor_fmt = dict(zip(format_fields, tumor_values))
            t_a, t_b = parse_copy_number(tumor_fmt.get(cn_field, ".,."))

            # Require numeric copy numbers for all alleles; otherwise skip variant
            if not (n_a.isnumeric() and n_b.isnumeric() and t_a.isnumeric() and t_b.isnumeric()):
                continue

            written_variants += 1
            rows.append(
                {
                    "sample": sample_name,
                    "tmp_id": tmp_id,
                    "tumor_allele_a": t_a,
                    "tumor_allele_b": t_b,
                    "normal_allele_a": n_a,
                    "normal_allele_b": n_b,
                    "end": end,
                    "period": period,
                    "ref": ref,
                    "motif": motif,
                    "genotype_separator": genotype_separator,
                }
            )

    # Simple logging
    if filter_by_pass:
        denom = written_variants + filter_not_passed_count
        if denom > 0:
            pct = filter_not_passed_count / denom * 100
            print(
                f"{path.name}: skipped {filter_not_passed_count} "
                f"({pct:.2f}%) variants due to FILTER != PASS"
            )

    if filter_by_perfect:
        total_checked = written_variants + non_perfect_count
        if total_checked > 0:
            pct = non_perfect_count / total_checked * 100
            print(
                f"{path.name}: skipped {non_perfect_count} "
                f"({pct:.2f}%) variants due to PERFECT == FALSE"
            )

    return rows


def parse_vcf_files(
    input_dir: str | Path,
    *,
    filter_by_pass: bool = True,
    filter_by_perfect: bool = True,
) -> pd.DataFrame:
    """
    Process all VCF (``.vcf`` / ``.vcf.gz``) files in a directory into a DataFrame.

    Supports GangSTR and conSTRain STR-annotated VCFs, as well as VCFs annotated
    with ``strvcf_annotator``.

    If a file causes an error, it is skipped and a message is printed.

    Parameters
    ----------
    input_dir : str or pathlib.Path
        Directory containing ``.vcf`` or ``.vcf.gz`` files.
    filter_by_pass : bool, optional
        If ``True``, keep only records with ``FILTER == "PASS"``.
    filter_by_perfect : bool, optional
        If ``True``, keep only records with ``INFO/PERFECT != "FALSE"`` when present.

    Returns
    -------
    pandas.DataFrame
        Parsed STR records concatenated across all input files.

        Columns:

        - ``sample``
        - ``tmp_id``
        - ``tumor_allele_a``
        - ``tumor_allele_b``
        - ``normal_allele_a``
        - ``normal_allele_b``
        - ``end``
        - ``period``
        - ``ref``
        - ``motif``
        - ``genotype_separator``
    """
    input_dir = Path(input_dir)
    all_rows: list[dict] = []

    for file in sorted(input_dir.iterdir()):
        # accept foo.vcf and foo.vcf.gz
        if file.name.endswith(".vcf") or file.name.endswith(".vcf.gz"):
            print(f"Processing {file.name}...")
            try:
                rows = process_vcf_to_rows(
                    file,
                    filter_by_pass=filter_by_pass,
                    filter_by_perfect=filter_by_perfect,
                )
            except Exception as e:
                # skip problematic file but continue with the rest
                print(f"Skipping {file.name} due to error: {e}")
                continue

            all_rows.extend(rows)

    columns = [
        "sample",
        "tmp_id",
        "tumor_allele_a",
        "tumor_allele_b",
        "normal_allele_a",
        "normal_allele_b",
        "end",
        "period",
        "ref",
        "motif",
        "genotype_separator",
    ]

    if not all_rows:
        return pd.DataFrame(columns=columns)

    return pd.DataFrame(all_rows, columns=columns)


def save_counts_matrix(mutations_data: pd.DataFrame, output_csv: str | Path):
    """
    Save a mutation counts matrix to a CSV file.

    Parameters
    ----------
    mutations_data : pandas.DataFrame
        DataFrame containing mutation count data to be written to disk.
    output_csv : str or pathlib.Path
        Path to the output CSV file.

    Returns
    -------
    None
    """
    mutations_data.to_csv(output_csv, index=False)
