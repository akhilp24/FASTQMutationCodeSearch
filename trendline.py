import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

def plot_trendlines(data, output_path):
    # Set up the plotting style
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    
    # Create a 2x2 subplot layout
    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    fig.suptitle('Mutation Rates vs Age', fontsize=16, fontweight='bold')
    
    # Define the variables to plot
    variables = [
        'total_mutations_over_total_g_strand_2xrepeats_per_1k',
        'G_T_g1_per_1k',  # pos1 mutation rate
        'G_T_g2_per_1k',  # pos2 mutation rate  
        'G_T_g3_per_1k'   # pos3 mutation rate
    ]
    
    titles = [
        'Total Mutations per 1000bp',
        'G < T at Position 1 Mutation Rate per 1000bp',
        'G < T at Position 2 Mutation Rate per 1000bp',
        'G < T at Position 3 Mutation Rate per 1000bp'
    ]
    
    # Plot each variable
    for i, (var, title) in enumerate(zip(variables, titles)):
        row = i // 2
        col = i % 2
        ax = axes[row, col]
        
        # Remove rows with missing Age values
        plot_data = data.dropna(subset=['Age', var])
        
        if len(plot_data) > 0:
            # Create scatter plot
            sns.scatterplot(data=plot_data, x='Age', y=var, ax=ax, alpha=0.6)
            
            # Add trendline
            if len(plot_data) > 1:
                sns.regplot(data=plot_data, x='Age', y=var, ax=ax, 
                          scatter=False, line_kws={'color': 'blue', 'linestyle': '--'})
            
            ax.set_title(title, fontweight='bold')
            ax.set_xlabel('Age')
            ax.set_ylabel('Mutations per 1000bp')
        else:
            ax.text(0.5, 0.5, 'No data available', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=12)
            ax.set_title(title, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def plot_spearman_correlations(data, output_path):
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    
    variables = [
        'total_mutations_over_total_g_strand_2xrepeats_per_1k',
        'G_T_g1_per_1k',
        'G_T_g2_per_1k',
        'G_T_g3_per_1k'
    ]
    titles = [
        'Total Mutations per 1000bp',
        'G < T at Position 1 Mutation Rate per 1000bp',
        'G < T at Position 2 Mutation Rate per 1000bp',
        'G < T at Position 3 Mutation Rate per 1000bp'
    ]
    
    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    fig.suptitle("Spearman's Rank Correlation: Mutation Rates vs Age", fontsize=16, fontweight='bold')
    
    for i, (var, title) in enumerate(zip(variables, titles)):
        row = i // 2
        col = i % 2
        ax = axes[row, col]
        plot_data = data.dropna(subset=['Age', var])
        if len(plot_data) > 1:
            corr, pval = spearmanr(plot_data['Age'], plot_data[var])
            sns.scatterplot(data=plot_data, x='Age', y=var, ax=ax, alpha=0.6)
            ax.set_title(f"{title}\nSpearman r = {corr:.2f}, p = {pval:.2g}", fontweight='bold')
            ax.set_xlabel('Age')
            ax.set_ylabel('Mutations per 1000bp')
        else:
            ax.text(0.5, 0.5, 'No data available', ha='center', va='center', 
                    transform=ax.transAxes, fontsize=12)
            ax.set_title(title, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def plot_trendlines_main():
    data = pd.read_csv("telomere_analysis.csv")
    # Plot the trendlines
    plot_trendlines(data, "trendline.png")
    # Plot the Spearman correlation graph
    plot_spearman_correlations(data, "spearman_correlation.png")
    print("Trendline plot saved as 'trendline.png'")
    print("Spearman correlation plot saved as 'spearman_correlation.png'")

if __name__ == "__main__":
    plot_trendlines_main()