import os
import argparse
import pandas as pd
from SigProfilerExtractor import sigpro as sig

def load_matrix(path):
    for sep in [",", "\t"]:
        try:
            df = pd.read_csv(path, index_col=0, sep=sep)
            if df.select_dtypes(include='number').shape[1] == df.shape[1]:
                return df, sep
        except Exception:
            continue
    raise ValueError(f"Unable to read {path} as a valid numeric mutation matrix.")

def describe_matrix(df, path):
    print(f"\n Matrix check passed for: {path}")
    print(f"Number of mutation types (signatures): {df.shape[0]}")
    print(f"Number of samples: {df.shape[1]}")
    print(f"Example mutation types: {list(df.index[:5])}")
    print(f"Example samples: {list(df.columns[:5])}\n")

def run_sigprofiler(input_files, output_dir, use_gpu, cpu, reference_genome):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for input_file in input_files:
        if not os.path.exists(input_file):
            print(f"⚠️  Warning: {input_file} not found. Skipping...")
            continue

        try:
            df, sep_used = load_matrix(input_file)
        except ValueError as e:
            print(f"Error: {e}")
            continue

        describe_matrix(df, input_file)

        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_path = os.path.join(output_dir, base_name)

        print(f"Running SigProfilerExtractor on {input_file} → Output: {output_path} (GPU: {use_gpu})")

        sig.sigProfilerExtractor(
            input_type="matrix",
            reference_genome=reference_genome,
            output=output_path,
            input_data=df,
            collapse_to_SBS96=False,
            make_decomposition_plots=True,
            context_type="custom",
            cpu=cpu if not use_gpu else 1,
            gpu=use_gpu
        )

    print("NMF Signature Extraction completed.")

def collect_input_files(input_path):
    if os.path.isfile(input_path):
        return [input_path]
    elif os.path.isdir(input_path):
        return [
            os.path.join(input_path, f)
            for f in os.listdir(input_path)
            if f.endswith(".csv") or f.endswith(".tsv")
        ]
    else:
        raise FileNotFoundError(f"{input_path} does not exist or is not valid.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SigProfilerExtractor on one or more STR count matrices.")
    parser.add_argument("--input_path", required=True, help="Path to input .csv/.tsv file or folder containing files")
    parser.add_argument("--output_dir", required=True, help="Directory where output will be saved")
    parser.add_argument("--gpu", action="store_true", help="Use GPU for signature extraction")
    parser.add_argument("--cpu", type=int, default=15, help="Number of CPUs to use (default: 15)")
    parser.add_argument("--reference_genome", default="GRCh38", help="Reference genome (default: GRCh38)")

    args = parser.parse_args()
    input_files = collect_input_files(args.input_path)
    run_sigprofiler(input_files, args.output_dir, use_gpu=args.gpu, cpu=args.cpu, reference_genome=args.reference_genome)
