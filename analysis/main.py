import os
import sys

# Patterns file used by generate_csv, plotting, and trendline. Set before other steps.
from patterns_config import set_patterns_file
_patterns_path = os.path.join(os.path.dirname(__file__), 'telomere_patterns_2x.json')
set_patterns_file(_patterns_path)

# Step 1: Generate CSV
print("[1/4] Generating telomere_analysis.csv ...")
from generate_csv import generate_csv
# data_dir = "/Users/akhilpeddikuppa/FieldLab/GreiderCodeSearch/greider_data_download"
# data_dir = "/Users/akhilpeddikuppa/FieldLab/GreiderCodeSearch/greider_data_download"
# generate_csv(data_dir)

# Step 2: Plot histograms
# print("[2/4] Plotting histograms ...")
# from histogram import plot_histograms_main
# plot_histograms_main()

# Step 3: Plot mutational signatures
print("[3/4] Plotting mutational signatures ...")
from plotting import plot_mutational_signatures_main
from plotting import plot_spearman_with_age_main
from plotting import plot_composite_score_main
# from plotting import plot_mutation_r_heatmap_main
from plotting import plot_pairwise_r_heatmap_main
# plot_mutational_signatures_main()
# plot_spearman_with_age_main()
# plot_composite_score_main()
# plot_mutation_r_heatmap_main()
# plot_pairwise_r_heatmap_main()

# Step 4: Plot trendlines and correlations (all config defined here)
print("[4/4] Plotting trendlines and correlations ...")
_analysis_dir = os.path.dirname(__file__)
_TRENDLINE_INPUT_CSV = os.path.join(_analysis_dir, "telomere_analysis_2x_repeat_feb4_2026.csv")
_TRENDLINE_OUTPUT = os.path.join(_analysis_dir, "trendline.png")
_SPEARMAN_OUTPUT = os.path.join(_analysis_dir, "spearman_correlation.png")
_TRENDLINE_VARIABLES = [
    "total_mutations_over_total_g_strand_2xrepeats_per_1k",
    "g_strand_T>C_sum_per_1k",
    "g_strand_G>T_sum_per_1k",
    "g_strand_T>G_sum_per_1k",
]
_TRENDLINE_TITLES = [
    "Total Mutations per 1000bp",
    "T > C Mutation Rate per 1000bp",
    "G > T Mutation Rate per 1000bp",
    "T > G Mutation Rate per 1000bp",
]
from trendline import plot_trendlines_main
plot_trendlines_main(
    csv_path=_TRENDLINE_INPUT_CSV,
    trendline_output_path=_TRENDLINE_OUTPUT,
    spearman_output_path=_SPEARMAN_OUTPUT,
    variables=_TRENDLINE_VARIABLES,
    titles=_TRENDLINE_TITLES,
)

print("\nAll processing complete.")
