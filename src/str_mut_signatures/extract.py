import os
import gzip
import pandas as pd


def parse_info(info_field: str) -> dict:
    info = {}
    for item in info_field.split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            info[key] = val
    return info


def parse_repcn(repcn_str: str):
    parts = repcn_str.split(',')
    if len(parts) == 2:
        return parts[0], parts[1]
    elif len(parts) == 1:
        return parts[0], parts[0]
    else:
        return '.', '.'


def _open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def process_vcf_to_rows(path: str):
    """
    Parse a single VCF and return a list of dict rows.

    Assumes:
    - At least two samples in the VCF.
    - FIRST sample = NORMAL, SECOND sample = TUMOR.
    - STR annotations present in INFO/FORMAT.
    """
    rows = []
    filter_not_passed_count = 0
    written_variants = 0

    sample_name = os.path.basename(path).replace(".vcf", "").replace(".vcf.gz", "")
    header_cols = None
    normal_idx = None
    tumor_idx = None
    format_fields = None

    with _open_maybe_gzip(path) as f:
        for line in f:
            if line.startswith('#CHROM'):
                header_cols = line.strip().split('\t')
                sample_cols = header_cols[9:]
                if len(sample_cols) < 2:
                    raise ValueError(f"Expected at least 2 samples (normal, tumor) in VCF: {sample_cols}")
                # By convention: first = NORMAL, second = TUMOR
                normal_idx = 9
                tumor_idx = 10
                continue

            if line.startswith('#'):
                continue

            cols = line.strip().split('\t')
            if len(cols) < 10:
                continue  # malformed

            # FILTER
            if cols[6] != "PASS":
                filter_not_passed_count += 1
                continue

            info = parse_info(cols[7])

            # Keep only perfect calls if annotated
            if info.get("PERFECT", "") == "FALSE":
                continue

            chrom = cols[0]
            pos = cols[1]
            tmp_id = f"{chrom}_{pos}"

            end = info.get('END', '')
            period = info.get('PERIOD', '')
            ref = info.get('REF', '')
            motif = info.get('RU', '')

            # FORMAT & samples
            if format_fields is None:
                format_fields = cols[8].split(':')

            def get_repcn(sample_str: str):
                values = sample_str.split(':')
                fmt = dict(zip(format_fields, values))
                return parse_repcn(fmt.get('REPCN', '.,.'))

            try:
                normal_sample = cols[normal_idx]
                tumor_sample = cols[tumor_idx]
            except IndexError:
                raise ValueError(f"VCF {path} does not contain expected normal/tumor sample columns.")

            n_a, n_b = get_repcn(normal_sample)
            t_a, t_b = get_repcn(tumor_sample)

            # Require numeric REPCN for all alleles
            if not (n_a.isnumeric() and n_b.isnumeric() and t_a.isnumeric() and t_b.isnumeric()):
                continue

            written_variants += 1
            rows.append(
                {
                    'sample': sample_name,
                    'tmp_id': tmp_id,
                    'tumor_allele_a': t_a,
                    'tumor_allele_b': t_b,
                    'normal_allele_a': n_a,
                    'normal_allele_b': n_b,
                    'end': end,
                    'period': period,
                    'ref': ref,
                    'motif': motif,
                }
            )

    denom = written_variants + filter_not_passed_count
    if denom > 0:
        filter_not_passed_pt = filter_not_passed_count / denom * 100
        print(
            f"{os.path.basename(path)}: skipped {filter_not_passed_count} "
            f"({filter_not_passed_pt:.2f}%) variants due to FILTER != PASS"
        )

    return rows


def parse_vcf_files(input_dir: str) -> pd.DataFrame:
    """
    Process all VCF(.gz) files in a directory into a single DataFrame.

    Returns
    -------
    pandas.DataFrame
        Columns:
        ['sample', 'tmp_id',
         'tumor_allele_a', 'tumor_allele_b',
         'normal_allele_a', 'normal_allele_b',
         'end', 'period', 'ref', 'motif']
    """
    all_rows = []
    for file in os.listdir(input_dir):
        if file.endswith(".vcf") or file.endswith(".vcf.gz"):
            full_path = os.path.join(input_dir, file)
            print(f"Processing {file}...")
            try:
                rows = process_vcf_to_rows(full_path)
                all_rows.extend(rows)
            except Exception as e:
                print(f"Skipping {file} due to error: {e}")

    columns = [
        'sample', 'tmp_id',
        'tumor_allele_a', 'tumor_allele_b',
        'normal_allele_a', 'normal_allele_b',
        'end', 'period', 'ref', 'motif',
    ]

    if not all_rows:
        return pd.DataFrame(columns=columns)

    return pd.DataFrame(all_rows, columns=columns)


def save_counts_matrix(mutations_data, output_csv):
    mutations_data.to_csv(output_csv)
