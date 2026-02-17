import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Read the CSV file
csv_path = 'telomere_analysis_2x_repeat.csv'
df = pd.read_csv(csv_path)

# Clean and filter data
# Some columns may have missing or non-numeric values, so we coerce errors
for col in ['Age', 'Telomere_Length', 'total_mutations_over_total_g_strand_2xrepeats_per_1k']:
    if col not in df.columns:
        raise ValueError(f"Column '{col}' not found in CSV.")
    df[col] = pd.to_numeric(df[col], errors='coerce')

df = df.dropna(subset=['Age', 'Telomere_Length', 'total_mutations_over_total_g_strand_2xrepeats_per_1k'])

# Prepare data for plotting
ages = df['Age']
telomere_lengths = df['Telomere_Length']
mutations = df['total_mutations_over_total_g_strand_2xrepeats_per_1k']

plt.figure(figsize=(10, 7))
# Create a custom colormap from green (low mutations) to red (high mutations)
from matplotlib.colors import LinearSegmentedColormap
colors = ['green', 'yellow', 'red']
n_bins = 100
cmap = LinearSegmentedColormap.from_list('green_to_red', colors, N=n_bins)
sc = plt.scatter(ages, telomere_lengths, c=mutations, cmap=cmap, s=60, edgecolor='k', alpha=0.8)
cbar = plt.colorbar(sc)
cbar.set_label('total_mutations_over_total_g_strand_2xrepeats_per_1k')

plt.xlabel('Age')
plt.ylabel('Telomere Length')
plt.title('Telomere Length vs Age (colored by total_mutations_over_total_g_strand_2xrepeats_per_1k)')

# Highlight and label top 3 outliers (highest mutation rates)
outlier_indices = np.argsort(mutations)[-3:]
for idx in outlier_indices:
    plt.annotate(
        df.iloc[idx]['FileName'] if 'FileName' in df.columns else str(idx),
        (ages.iloc[idx], telomere_lengths.iloc[idx]),
        textcoords="offset points", xytext=(5,5), ha='left', color='red', fontsize=9, weight='bold',
        arrowprops=dict(arrowstyle='->', color='red', lw=1)
    )

plt.tight_layout()
plt.savefig('scatterplot.png')
