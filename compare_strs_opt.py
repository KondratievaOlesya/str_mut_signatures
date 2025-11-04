import os
import csv

def parse_info(info_field):
    info = {}
    for item in info_field.split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            info[key] = val
    return info

def parse_repcn(repcn_str):
    parts = repcn_str.split(',')
    if len(parts) == 2:
        return parts[0], parts[1]
    elif len(parts) == 1:
        return parts[0], parts[0]
    else:
        return '.', '.'

def process_vcf(path, writer):
    sample_name = os.path.basename(path).replace(".vcf", "").replace(".vcf.gz", "")
    tumor_sample_id = os.path.basename(path).split('.')[2]  # Get SAxxxxxx
    filter_not_passed_count = 0
    written_variants = 0
    with open(path) as f:
        for line in f:
            if line.startswith('#CHROM'):
                header_cols = line.strip().split('\t')
                sample_cols = header_cols[9:]
                if len(sample_cols) < 2:
                    raise ValueError(f"Expected at least 2 samples in VCF: {sample_cols}")

                try:
                    tumor_idx = sample_cols.index(tumor_sample_id) + 9
                    normal_idx = 9 if tumor_idx == 10 else 10
                except ValueError:
                    raise ValueError(f"Tumor sample {tumor_sample_id} not found in VCF header: {sample_cols}")
                format_fields = None
                continue

            if line.startswith('#'):
                continue

            cols = line.strip().split('\t')
            if len(cols) < 10:
                continue  # malformed line
            filter_col = cols[6]
            if filter_col != "PASS":
                filter_not_passed_count += 1
                continue # Skip low quality variants
            chrom = cols[0]
            pos = cols[1]
            tmp_id = f"{chrom}_{pos}"

            info = parse_info(cols[7])
            if info.get("PERFECT", '') == "FALSE":
                continue
            end = info.get('END', '')
            period = info.get('PERIOD', '')
            ref = info.get('REF', '')
            motif = info.get('RU', '')

            if format_fields is None:
                format_fields = cols[8].split(':')

            tumor_sample = cols[tumor_idx]
            normal_sample = cols[normal_idx]

            def get_repcn(sample_str):
                values = sample_str.split(':')
                fmt = dict(zip(format_fields, values))
                return parse_repcn(fmt.get('REPCN', '.,.'))

            t_a, t_b = get_repcn(tumor_sample)
            n_a, n_b = get_repcn(normal_sample)
            if not (t_a.isnumeric() and t_b.isnumeric() and n_a .isnumeric() and n_b.isnumeric()):
                continue
            written_variants += 1
            row = {
                'sample': sample_name,
                'tmp_id': tmp_id,
                'tumor_allele_a': t_a,
                'tumor_allele_b': t_b,
                'normal_allele_a': n_a,
                'normal_allele_b': n_b,
                'end': end,
                'period': period,
                'ref': ref,
                'motif': motif
            }
            writer.writerow(row)
    filter_not_passed_pt = filter_not_passed_count / (written_variants + filter_not_passed_count) * 100
    print(f"Skip {filter_not_passed_count} ({filter_not_passed_pt}%) variants due to filter not PASS")


def main():
    input_dir = "/media/pho/Elements1/muto_intl_str"
    output_csv = "/home/pho/str_signatures/data/str_df/muto_intl_mutect2_pass.csv"

    fieldnames = ['sample', 'tmp_id', 'tumor_allele_a', 'tumor_allele_b',
                  'normal_allele_a', 'normal_allele_b', 'end', 'period', 'ref', 'motif']

    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for file in os.listdir(input_dir):
            if file.endswith(".vcf"):
                full_path = os.path.join(input_dir, file)
                print(f"Processing {file}...")
                try:
                    process_vcf(full_path, writer)
                except Exception as e:
                    print(f"Skipping {file} due to error: {e}")

    print(f"Saved to {output_csv}")


if __name__ == "__main__":
    main()
