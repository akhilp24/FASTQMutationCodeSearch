import argparse
import os
from datetime import datetime

from generate_csv import generate_csv
from plotting import (
    plot_mutational_signatures_main,
    plot_spearman_with_age_main,
    plot_pairwise_r_heatmap_main,
    plot_trendlines_main,
    plot_histograms_from_csv,
    _default_csv_path_from_patterns,
)

_ANALYSIS_DIR = os.path.dirname(__file__)
_TRENDLINES_DIR = os.path.join(_ANALYSIS_DIR, 'trendlines')
_SPEARMAN_CORRELATIONS_DIR = os.path.join(_ANALYSIS_DIR, 'spearman_correlations')
_HISTOGRAM_DIR = os.path.join(_ANALYSIS_DIR, 'histograms')


def _unique_output_path(directory: str, base_name: str, ext: str) -> str:
    """Create directory if needed; return path with unique filename (timestamp) to avoid duplicates."""
    os.makedirs(directory, exist_ok=True)
    timestamp = datetime.now().strftime('%Y-%m-%d_%H%M%S')
    filename = f"{base_name}_{timestamp}.{ext}"
    return os.path.join(directory, filename)


def _parse_args():
    parser = argparse.ArgumentParser(
        description="CLI for generating telomere mutation CSVs and running analysis plots."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # generate: build telomere_analysis_*.csv from FASTQ/FASTA + patterns + metadata
    gen = subparsers.add_parser(
        "generate",
        help="Generate telomere_analysis CSV from FASTQ/FASTA files.",
    )
    gen.add_argument(
        "--patterns",
        required=True,
        help="Path to telomere patterns JSON (e.g. telomere_patterns_2x.json).",
    )
    gen.add_argument(
        "--metadata",
        required=False,
        help="Optional path to greider_methods_table_s2_outliers_removed.csv (age/length table).",
    )
    gen.add_argument(
        "--fastq-dir",
        required=True,
        help="Directory containing FASTQ/FASTA files.",
    )
    gen.add_argument(
        "--csv-out",
        help=(
            "Optional explicit output CSV path. "
            "If omitted, generate_csv will derive a filename from the patterns version."
        ),
    )

    # plot: run plotting on an existing CSV
    plot = subparsers.add_parser(
        "plot",
        help="Run plotting pipelines on an existing telomere_analysis CSV.",
    )
    plot.add_argument(
        "--patterns",
        required=True,
        help="Path to telomere patterns JSON used to generate the CSV (for version labels).",
    )
    plot.add_argument(
        "--csv",
        help=(
            "Path to telomere_analysis CSV. If omitted, "
            "it will be inferred from the patterns version."
        ),
    )
    plot.add_argument(
        "--no-hist",
        action="store_true",
        help="Skip histogram / per-file mutation plots.",
    )
    plot.add_argument(
        "--no-signatures",
        action="store_true",
        help="Skip mutational signature plots.",
    )
    plot.add_argument(
        "--no-spearman",
        action="store_true",
        help="Skip simple Spearman vs Age scatter plots.",
    )
    plot.add_argument(
        "--no-pairwise",
        action="store_true",
        help="Skip pairwise Spearman heatmap.",
    )
    plot.add_argument(
        "--no-trendlines",
        action="store_true",
        help="Skip 2x2 trendline and Spearman correlation grids.",
    )
    plot.add_argument(
        "--no-curve",
        action="store_true",
        help="Skip curve-fitting analysis.",
    )

    # run: convenience command that generates CSV then runs plotting
    run = subparsers.add_parser(
        "run",
        help="Generate CSV from FASTQ/FASTA and then run plotting in one command.",
    )
    run.add_argument(
        "--patterns",
        required=True,
        help="Path to telomere patterns JSON (e.g. telomere_patterns_2x.json).",
    )
    run.add_argument(
        "--metadata",
        required=False,
        help="Optional path to greider_methods_table_s2_outliers_removed.csv (age/length table).",
    )
    run.add_argument(
        "--fastq-dir",
        required=True,
        help="Directory containing FASTQ/FASTA files.",
    )
    run.add_argument(
        "--csv-out",
        help=(
            "Optional explicit output CSV path. "
            "If omitted, generate_csv will derive a filename from the patterns version."
        ),
    )
    run.add_argument(
        "--no-hist",
        action="store_true",
        help="Skip histogram / per-file mutation plots.",
    )
    run.add_argument(
        "--no-signatures",
        action="store_true",
        help="Skip mutational signature plots.",
    )
    run.add_argument(
        "--no-spearman",
        action="store_true",
        help="Skip simple Spearman vs Age scatter plots.",
    )
    run.add_argument(
        "--no-pairwise",
        action="store_true",
        help="Skip pairwise Spearman heatmap.",
    )
    run.add_argument(
        "--no-trendlines",
        action="store_true",
        help="Skip 2x2 trendline and Spearman correlation grids.",
    )
    run.add_argument(
        "--no-curve",
        action="store_true",
        help="Skip curve-fitting analysis.",
    )

    return parser.parse_args()


def _resolve_csv_path(patterns_path: str, explicit_csv: str | None) -> str:
    if explicit_csv:
        return explicit_csv
    # Use plotting helper to infer the filename from the patterns version
    return _default_csv_path_from_patterns(patterns_path)


def _run_generate(args) -> str:
    print("[1/2] Generating telomere_analysis CSV ...")
    csv_out = args.csv_out
    generate_csv(
        data_dir=args.fastq_dir,
        metadata_file_path=args.metadata,
        patterns_file_path=args.patterns,
        output_csv_path=csv_out,
    )
    # If the caller did not specify csv_out, reconstruct what generate_csv used.
    csv_path = csv_out or _default_csv_path_from_patterns(args.patterns)
    print(f"CSV written to: {csv_path}")
    return csv_path


def _run_plots(args, csv_path: str):
    print("[2/2] Running plots ...")

    if not args.no_hist:
        os.makedirs(_HISTOGRAM_DIR, exist_ok=True)
        plot_histograms_from_csv(csv_path, output_dir=_HISTOGRAM_DIR)

    if not args.no_signatures:
        plot_mutational_signatures_main(args.patterns)

    if not args.no_spearman:
        plot_spearman_with_age_main(args.patterns)

    if not args.no_pairwise:
        plot_pairwise_r_heatmap_main(args.patterns)

    if not args.no_trendlines:
        trendline_output = _unique_output_path(_TRENDLINES_DIR, "trendline", "png")
        spearman_output = _unique_output_path(_SPEARMAN_CORRELATIONS_DIR, "spearman_correlation", "png")
        _TRENDLINE_VARIABLES = [
            "total_mutations_over_total_g_strand_per_1k",
            "g_strand_T>C_sum_per_1k",
            "g_strand_G>T_sum_per_1k",
            "g_strand_T>G_sum_per_1k",
        ]
        _TRENDLINE_TITLES = [
            "Total Mutations Normalized",
            "T > C Mutation Rate Normalized",
            "G > T Mutation Rate Normalized",
            "T > G Mutation Rate Normalized",
        ]
        plot_trendlines_main(
            csv_path=csv_path,
            trendline_output_path=trendline_output,
            spearman_output_path=spearman_output,
            variables=_TRENDLINE_VARIABLES,
            titles=_TRENDLINE_TITLES,
            patterns_file_path=args.patterns,
        )

    if not args.no_curve:
        # curve_fitting_analysis_main reads its CSV internally via _default_csv_path_from_patterns,
        # but that will match csv_path as long as you followed the generate/run flow.
        from plotting import curve_fitting_analysis_main

        curve_fitting_analysis_main(args.patterns)

    print("\nAll plotting complete.")


def main():
    args = _parse_args()

    if args.command == "generate":
        _run_generate(args)
    elif args.command == "plot":
        csv_path = _resolve_csv_path(args.patterns, args.csv)
        _run_plots(args, csv_path)
    elif args.command == "run":
        csv_path = _run_generate(args)
        _run_plots(args, csv_path)
    else:
        raise SystemExit(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
