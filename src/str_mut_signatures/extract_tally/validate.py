from __future__ import annotations

import gzip
from pathlib import Path
from typing import Literal, NamedTuple


class VCFValidationResult(NamedTuple):
    has_str_annotations: bool
    has_paired_samples: bool
    normal_sample: str | None
    tumor_sample: str | None
    # "|"  -> phased / directed GT  (e.g. 1|0)
    # "/"  -> unphased GT          (e.g. 1/0)
    # None -> GT not found / cannot determine from first record
    genotype_separator: str | None
    # Name of the FORMAT field giving allele copy number:
    # - "REPCN" for GangSTR
    # - "REPLEN" for conSTRain
    # - None if not detected
    copy_number_field: str | None
    # Best guess of STR caller based on header fields:
    # - "gangstr"
    # - "constrain"
    # - "unknown" if we can't decide
    caller: Literal["gangstr", "constrain", "unknown"]


# Required header fields for STR-annotated VCFs (common to all)
REQUIRED_INFO_FIELDS = {"RU", "REF"}

# Possible FORMAT fields that encode allele-level repeat copy number
COPY_NUMBER_FIELDS = {"REPCN", "REPLEN"}


def open_maybe_gzip(path: Path):
    """
    Open .vcf or .vcf.gz transparently in text mode.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def validate_vcf(vcf_path: str | Path) -> VCFValidationResult:
    """
    Inspect only the VCF header to determine:

    - whether required STR annotations are defined in header for
      GangSTR or conSTRain:
        * INFO/RU
        * INFO/REF
        * FORMAT/REPCN (GangSTR) or FORMAT/REPLEN (conSTRain)
    - whether there are at least two samples (normal, tumor)
      in the #CHROM header line.
      By convention:
        * first sample column after FORMAT = normal
        * second sample column after FORMAT = tumor
    - which copy-number field is present (REPCN or REPLEN)
    - best guess of the STR caller (gangstr / constrain / unknown)

    Parameters
    ----------
    vcf_path:
        Path to a .vcf or .vcf.gz file.

    Returns
    -------
    VCFValidationResult
        has_str_annotations : True if RU/REF are present and at least one
                              of REPCN / REPLEN exists in FORMAT.
        has_paired_samples  : True if at least two sample columns present.
        normal_sample       : Name of the first sample (or None).
        tumor_sample        : Name of the second sample (or None).
        genotype_separator  : "|" for phased, "/" for unphased, or None if
                              GT could not be determined from the first record.
        copy_number_field   : "REPCN", "REPLEN", or None.
        caller              : "gangstr", "constrain", or "unknown".
    """
    path = Path(vcf_path)

    info_ids: set[str] = set()
    format_ids: set[str] = set()
    sample_names: list[str] = []
    format_col_idx: int | None = None
    first_record_line: str | None = None

    with open_maybe_gzip(path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")

            if line.startswith("##INFO=<ID="):
                # Example: ##INFO=<ID=RU,Number=1,Type=String,...>
                try:
                    after_id = line.split("ID=", 1)[1]
                    id_val = after_id.split(",", 1)[0]
                    info_ids.add(id_val)
                except IndexError:
                    continue

            elif line.startswith("##FORMAT=<ID="):
                # Example: ##FORMAT=<ID=REPCN,Number=R,Type=Integer,...>
                try:
                    after_id = line.split("ID=", 1)[1]
                    id_val = after_id.split(",", 1)[0]
                    format_ids.add(id_val)
                except IndexError:
                    continue

            elif line.startswith("#CHROM"):
                # Standard VCF header columns:
                # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE...]
                fields = line.lstrip("#").split()
                try:
                    format_col_idx = fields.index("FORMAT")
                except ValueError:
                    format_col_idx = None
                    sample_names = []
                else:
                    sample_names = fields[format_col_idx + 1 :]

            elif line.startswith("#"):
                # other header lines, ignore
                continue

            else:
                # first non-header line: first variant record
                first_record_line = line
                break

    # --- STR annotations / copy-number field detection -----------------------
    has_ru_ref = REQUIRED_INFO_FIELDS.issubset(info_ids)
    available_cn_fields = COPY_NUMBER_FIELDS.intersection(format_ids)

    if has_ru_ref and available_cn_fields:
        has_str_annotations = True
        # Prefer REPCN if both are present (very unlikely, but defined)
        copy_number_field = "REPCN" if "REPCN" in available_cn_fields else "REPLEN"
    else:
        has_str_annotations = False
        copy_number_field = None

    # --- Caller heuristic (GangSTR vs conSTRain) -----------------------------
    # GangSTR-specific clues:
    #  - FORMAT: REPCN, REPCI, RC, ENCLREADS, FLNKREADS, ML, INS, STDERR, QEXP, GGL
    #  - INFO: GRID, EXPTHRESH, STUTTERUP, STUTTERDOWN, STUTTERP
    gangstr_markers_info = {"GRID", "EXPTHRESH", "STUTTERUP", "STUTTERDOWN", "STUTTERP"}
    gangstr_markers_format = {
        "REPCI",
        "RC",
        "ENCLREADS",
        "FLNKREADS",
        "ML",
        "INS",
        "STDERR",
        "QEXP",
        "GGL",
    }

    # conSTRain-specific clues:
    #  - FORMAT: CN, FREQS, REPLEN
    constrain_markers_format = {"CN", "FREQS", "REPLEN"}

    is_gangstr = bool(
        ("REPCN" in format_ids)
        and (gangstr_markers_info & info_ids or gangstr_markers_format & format_ids)
    )
    is_constrain = bool(("REPLEN" in format_ids) and (constrain_markers_format & format_ids))

    if is_gangstr and not is_constrain:
        caller: Literal["gangstr", "constrain", "unknown"] = "gangstr"
    elif is_constrain and not is_gangstr:
        caller = "constrain"
    elif is_gangstr and is_constrain:
        # Extremely unlikely, but in case of mixed headers fall back to unknown
        caller = "unknown"
    else:
        caller = "unknown"

    # --- Sample info ---------------------------------------------------------
    has_paired_samples = len(sample_names) >= 2

    normal_sample = sample_names[0] if has_paired_samples else None
    tumor_sample = sample_names[1] if has_paired_samples else None

    # --- Detect GT separator ("|" vs "/") from the first record, if possible --
    genotype_separator: str | None = None

    if first_record_line is not None and format_col_idx is not None and has_paired_samples:
        fields = first_record_line.split()
        # Sanity: we need FORMAT + at least two samples
        if len(fields) > format_col_idx + 2:
            format_field = fields[format_col_idx]  # e.g. "GT:AD:DP:REPCN"
            format_keys = format_field.split(":")

            # try to find GT index within the FORMAT string
            try:
                gt_idx = format_keys.index("GT")
            except ValueError:
                gt_idx = None

            if gt_idx is not None:
                # sample columns: everything after FORMAT
                sample_fields = fields[format_col_idx + 1 :]

                # look at first two samples (normal, tumor)
                gt_tokens: list[str] = []
                for sf in sample_fields[:2]:
                    vals = sf.split(":")
                    if len(vals) > gt_idx:
                        gt_tokens.append(vals[gt_idx])

                # Decide separator based on any GT token
                for gt in gt_tokens:
                    if "|" in gt:
                        genotype_separator = "|"
                        break
                    if "/" in gt:
                        genotype_separator = "/"
                        # don't break yet in case we later see "|",
                        # but realistically they will be consistent

    return VCFValidationResult(
        has_str_annotations=has_str_annotations,
        has_paired_samples=has_paired_samples,
        normal_sample=normal_sample,
        tumor_sample=tumor_sample,
        genotype_separator=genotype_separator,
        copy_number_field=copy_number_field,
        caller=caller,
    )
