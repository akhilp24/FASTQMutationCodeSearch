import sys

# Step 1: Generate CSV
print("[1/4] Generating telomere_analysis.csv ...")
from generate_csv import generate_csv
generate_csv()

# Step 2: Plot histograms
print("[2/4] Plotting histograms ...")
from histogram import plot_histograms_main
plot_histograms_main()

# Step 3: Plot mutational signatures
print("[3/4] Plotting mutational signatures ...")
from plotting import plot_mutational_signatures_main
plot_mutational_signatures_main()

# Step 4: Plot trendlines and correlations
print("[4/4] Plotting trendlines and correlations ...")
from trendline import plot_trendlines_main
plot_trendlines_main()

print("\nAll processing complete.")
