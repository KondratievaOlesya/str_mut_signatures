import click
import os
import sys
from .extract import parse_vcf_files, save_counts_matrix
from .matrix_builder import build_mutation_matrix, create_comprehensive_matrix
from .nmf import run_nmf_decomposition, run_sigprofiler_extractor


@click.group()
def cli():
    """STR mutation signature analysis tool."""
    pass


@cli.command()
@click.option('--vcf-dir', required=True, help='Directory containing VCF files')
@click.option('--mode', default='ru',
              type=click.Choice(['ru', 'len', 'change_only', 'binned', 'all']),
              help='Matrix building mode')
@click.option('--out-matrix', required=True, help='Output counts matrix file')
@click.option('--out-dir', default='.', help='Output directory for multiple matrices')
def extract(vcf_dir, mode, out_matrix, out_dir):
    """Extract STR mutations from VCF files and create counts matrix."""

    if not os.path.exists(vcf_dir):
        raise click.ClickException(f"VCF directory does not exist: {vcf_dir}")

    click.echo(f"Parsing VCF files from: {vcf_dir}")
    click.echo(f"Using mode: {mode}")

    # Parse VCF files
    mutations_data = parse_vcf_files(vcf_dir)

    if not mutations_data:
        click.echo("No valid STR mutations found in the VCF files.")
        return

    click.echo(f"Found mutations for {len(mutations_data)} samples")

    # Build matrix/matrices
    if mode == 'all':
        matrices = create_comprehensive_matrix(mutations_data)

        os.makedirs(out_dir, exist_ok=True)
        for matrix_mode, df in matrices.items():
            if not df.empty:
                output_file = os.path.join(out_dir, f"counts_{matrix_mode}.tsv")
                save_counts_matrix(df, output_file)
                click.echo(f"Saved {matrix_mode} matrix: {output_file} (shape: {df.shape})")
    else:
        counts_df = build_mutation_matrix(mutations_data, mode=mode)

        if counts_df.empty:
            click.echo("No mutations to count after filtering.")
            return

        save_counts_matrix(counts_df, out_matrix)
        click.echo(f"Counts matrix saved to: {out_matrix}")
        click.echo(f"Matrix shape: {counts_df.shape}")
        click.echo(f"Mutation types: {len(counts_df.columns)}")


@cli.command()
@click.option('--matrix', required=True, help='Input counts matrix file')
@click.option('--outdir', required=True, help='Output directory for NMF results')
@click.option('--components', default=3, help='Number of components for NMF')
@click.option('--use-sigprofiler', is_flag=True, help='Use SigProfilerExtractor instead of sklearn')
def nmf(matrix, outdir, components, use_sigprofiler):
    """Run NMF decomposition on counts matrix."""

    if not os.path.exists(matrix):
        raise click.ClickException(f"Counts matrix file does not exist: {matrix}")

    click.echo(f"Running NMF decomposition on: {matrix}")
    click.echo(f"Output directory: {outdir}")
    click.echo(f"Number of components: {components}")

    if use_sigprofiler:
        click.echo("Using SigProfilerExtractor...")
        run_sigprofiler_extractor(matrix, outdir)
    else:
        click.echo("Using sklearn NMF...")
        run_nmf_decomposition(matrix, outdir, n_components=components)

    click.echo("NMF analysis completed successfully.")


if __name__ == '__main__':
    cli()