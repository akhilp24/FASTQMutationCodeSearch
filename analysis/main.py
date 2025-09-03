import sys

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

# Step 4: Plot trendlines and correlations
print("[4/4] Plotting trendlines and correlations ...")
from trendline import plot_trendlines_main
plot_trendlines_main()

print("\nAll processing complete.")
