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
        'total_mutations_over_total_g_strand_2xrepeats_per_1k',
        'g_strand_A>G_sum_per_1k',  # pos1 mutation rate
        'g_strand_T>G_sum_per_1k',  # pos2 mutation rate  
        'g_strand_T>C_sum_per_1k'   # pos3 mutation rate
    ]
    
    titles = [
        'Total Mutations Normalized',
        'G > A Mutation Rate Normalized',
        'T > G Mutation Rate Normalized',
        'T > C Mutation Rate Normalized'
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

def plot_mutations_per_file(data, output_path):
    """
    Create a histogram showing total number of mutations for each file.
    
    Args:
        data: DataFrame with mutation data
        output_path: Path to save the histogram
    """
    # Set up the plotting style
    sns.set_style("whitegrid")
    plt.figure(figsize=(16, 10))
    
    # Calculate total mutations for each file by summing all raw mutation columns
    # Get all mutation columns (those that contain 'mutations' but not 'per_1k')
    mutation_columns = [col for col in data.columns 
                       if 'mutations' in col and 'per_1k' not in col]
    
    # Calculate total mutations per file
    data_with_totals = data.copy()
    data_with_totals['Total_Mutations'] = data[mutation_columns].sum(axis=1)
    
    # Remove rows with missing data
    plot_data = data_with_totals.dropna(subset=['FileName', 'Total_Mutations'])
    
    if len(plot_data) > 0:
        # Sort by total mutations for better visualization
        plot_data = plot_data.sort_values('Total_Mutations', ascending=True)
        
        # Create the bar plot
        plt.figure(figsize=(16, 10))
        bars = plt.bar(range(len(plot_data)), plot_data['Total_Mutations'], 
                      color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.5)
        
        # Customize the plot
        plt.title('Total Number of Mutations per File', fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('Files', fontsize=12, fontweight='bold')
        plt.ylabel('Total Number of Mutations', fontsize=12, fontweight='bold')
        
        # Set x-axis labels to file names (rotated for readability)
        plt.xticks(range(len(plot_data)), plot_data['FileName'], rotation=45, ha='right')
        
        # Add value labels on top of bars for better readability
        for i, bar in enumerate(bars):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + max(plot_data['Total_Mutations'])*0.01,
                    f'{int(height)}', ha='center', va='bottom', fontsize=8)
        
        # Add grid for better readability
        plt.grid(axis='y', alpha=0.3)
        
        # Add some statistics as text
        mean_mutations = plot_data['Total_Mutations'].mean()
        median_mutations = plot_data['Total_Mutations'].median()
        max_mutations = plot_data['Total_Mutations'].max()
        min_mutations = plot_data['Total_Mutations'].min()
        
        stats_text = f'Mean: {mean_mutations:.0f}\nMedian: {median_mutations:.0f}\nMin: {min_mutations:.0f}\nMax: {max_mutations:.0f}'
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
                fontsize=10, verticalalignment='top', 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
    else:
        plt.text(0.5, 0.5, 'No data available', ha='center', va='center', 
                transform=plt.gca().transAxes, fontsize=12)
        plt.title('Total Number of Mutations per File', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    return len(plot_data)

def plot_histograms_main():
    data = pd.read_csv("telomere_analysis_2x_repeat.csv")
    
    # Plot the age-based histograms
    plot_histograms(data, "histogram.png")
    print("Histogram plot saved as 'histogram.png'")
    
    # Plot the mutations per file histogram
    num_files = plot_mutations_per_file(data, "mutations_per_file_histogram.png")
    print(f"Mutations per file histogram saved as 'mutations_per_file_histogram.png'")
    print(f"Processed {num_files} files")

if __name__ == "__main__":
    plot_histograms_main()
