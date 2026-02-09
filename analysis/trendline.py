import json
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, linregress


def _get_patterns_version(patterns_file_path):
    """Read version string from patterns JSON (path passed from main.py)."""
    try:
        with open(patterns_file_path) as f:
            return json.load(f).get('version', 'unknown')
    except Exception:
        return 'unknown'


def plot_trendlines(data, output_path, variables, titles, version):
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    
    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    fig.suptitle(f'Mutation Rates vs Age [{version}]', fontsize=16, fontweight='bold')
    
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
            
            # Add trendline and calculate R-squared
            if len(plot_data) > 1:
                sns.regplot(data=plot_data, x='Age', y=var, ax=ax, 
                          scatter=False, line_kws={'color': 'blue', 'linestyle': '--'})
                
                # Calculate R-squared value
                slope, intercept, r_value, p_value, std_err = linregress(plot_data['Age'], plot_data[var])
                r_squared = r_value ** 2
                
                # Add R-squared to title
                ax.set_title(f"{title}\nR² = {r_squared:.3f}", fontweight='bold')
            else:
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

def plot_spearman_correlations(data, output_path, variables, titles, version):
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    
    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    fig.suptitle(f"Spearman's Rank Correlation: Mutation Rates vs Age [{version}]", fontsize=16, fontweight='bold')
    
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

def plot_trendlines_main(
    csv_path,
    trendline_output_path,
    spearman_output_path,
    variables,
    titles,
    patterns_file_path,
):
    version = _get_patterns_version(patterns_file_path)
    data = pd.read_csv(csv_path)
    plot_trendlines(data, trendline_output_path, variables, titles, version)
    plot_spearman_correlations(data, spearman_output_path, variables, titles, version)
    print(f"Trendline plot saved as '{trendline_output_path}'")
    print(f"Spearman correlation plot saved as '{spearman_output_path}'")


if __name__ == "__main__":
    raise SystemExit(
        "This module is intended to be run via analysis/main.py so all inputs/outputs\n"
        "and variable/title lists live in one place.\n"
        "Run: python analysis/main.py"
    )