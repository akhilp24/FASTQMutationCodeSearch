import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_histograms(data, output_path):
    # Set up the plotting style
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    
    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    fig.suptitle('Mutation Rates by Age Groups (10-year bins)', fontsize=16, fontweight='bold')
    
    variables = [
        'total_mutations_over_total_g_strand_3xrepeats_per_1k',
        'g_strand_mutations_G>T_g1_per_1k',  # pos1 mutation rate
        'g_strand_mutations_G>T_g2_per_1k',  # pos2 mutation rate  
        'g_strand_mutations_G>T_g3_per_1k'   # pos3 mutation rate
    ]
    
    titles = [
        'Total Mutations per 1000bp',
        'G > T at Position 1 Mutation Rate per 1000bp',
        'G > T at Position 2 Mutation Rate per 1000bp',
        'G > T at Position 3 Mutation Rate per 1000bp'
    ]
    
    for i, (var, title) in enumerate(zip(variables, titles)):
        row = i // 2
        col = i % 2
        ax = axes[row, col]
        
        # Remove rows with missing Age values
        plot_data = data.dropna(subset=['Age', var])
        
        if len(plot_data) > 0:
            # Create age bins of 10 years
            age_bins = np.arange(0, plot_data['Age'].max() + 10, 10)
            age_labels = [f'{int(bin_start)}-{int(bin_start+9)}' for bin_start in age_bins[:-1]]
            
            # Add age group column
            plot_data['Age_Group'] = pd.cut(plot_data['Age'], bins=age_bins, labels=age_labels, include_lowest=True)
            
            # Create box plot with age groups on x-axis and mutation rates on y-axis
            sns.boxplot(data=plot_data, x='Age_Group', y=var, ax=ax)
            
            ax.set_title(title, fontweight='bold')
            ax.set_xlabel('Age Group (years)')
            ax.set_ylabel('Mutations per 1000bp')
            
            # Rotate x-axis labels for better readability
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
            
        else:
            ax.text(0.5, 0.5, 'No data available', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=12)
            ax.set_title(title, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def plot_histograms_main():
    data = pd.read_csv("telomere_analysis.csv")
    # Plot the histograms
    plot_histograms(data, "histogram.png")
    print("Histogram plot saved as 'histogram.png'")

if __name__ == "__main__":
    plot_histograms_main()
