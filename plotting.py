import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns

def plot_mutational_signature_row(row, mutation_types, mutation_columns, output_path):
    # Set seaborn style for better-looking plots
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    
    # Aggregate counts for each mutation type, position, and strand context for a single row
    bar_heights = []
    bar_colors = []
    bar_labels = []
    # Get all columns for total calculation
    all_columns = []
    for mut_type, contexts in mutation_columns.items():
        for context, cols in contexts.items():
            all_columns.extend(cols)
    total_mutations = row[all_columns].sum()
    for mut_label, color in mutation_types:
        contexts = mutation_columns[mut_label]
        for context_name, cols in contexts.items():
            for i, col in enumerate(cols):
                value = row[col]
                percentage = (value / total_mutations) * 100 if total_mutations > 0 else 0
                bar_heights.append(percentage)
                bar_colors.append(color)
                bar_labels.append(f"{mut_label} {context_name} pos{i+1}")
    
    x = np.arange(len(bar_heights))
    
    # Create figure with seaborn styling
    fig, ax = plt.subplots(figsize=(16, 10))
    
    # Create the bar plot with seaborn styling
    bars = sns.barplot(x=x, y=bar_heights, palette=bar_colors, ax=ax, edgecolor='black', linewidth=0.5)
    
    # Customize the plot
    for i, label in enumerate(bar_labels):
        ax.text(i, -max(bar_heights)*0.02, label, ha='center', va='center', 
                color='black', fontsize=9, fontweight='normal', rotation=45)
    
    ax.set_xticks([])
    ax.set_yticks(np.linspace(0, max(bar_heights), 6))
    ax.set_xlim(-0.5, len(bar_heights) - 0.5)
    ax.set_ylim(-max(bar_heights)*0.1, max(bar_heights) + max(bar_heights)*0.2)
    
    # Remove spines for cleaner look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Add grid with seaborn styling
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    ax.set_axisbelow(True)
    
    # Set labels and title with improved styling
    ax.set_ylabel('Percentage of Single Base Modifications', fontsize=14, fontweight='bold')
    
    # Title with age and filename
    age = row['Age'] if 'Age' in row else 'N/A'
    filename = row['FileName'] if 'FileName' in row else 'sample'
    title = f"Mutational Signatures by Position and Strand Context\nFile: {filename} | Age: {age} years"
    ax.set_title(title, fontsize=18, fontweight='bold', pad=30)
    
    # Improve layout
    plt.tight_layout(rect=[0, 0.15, 1, 0.95])
    
    # Save with high quality
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def plot_mutational_signatures(csv_path):
    # Set global seaborn style
    sns.set_theme(style="whitegrid", font_scale=1.1)
    
    df = pd.read_csv(csv_path)
    mutation_types = [
        ('C>A', 'blue'),
        ('C>G', 'black'),
        ('C>T', 'red'),
        ('G>A', 'gray'),
        ('G>C', 'green'),
        ('G>T', 'pink'),
    ]
    mutation_columns = {
        'C>A': {
            'C-strand': ['c_strand_mutations_C>A_c1', 'c_strand_mutations_C>A_c2', 'c_strand_mutations_C>A_c3'],
            'G-strand': ['g_strand_mutations_G>T_g1', 'g_strand_mutations_G>T_g2', 'g_strand_mutations_G>T_g3']
        },
        'C>G': {
            'C-strand': ['c_strand_mutations_C>G_c1', 'c_strand_mutations_C>G_c2', 'c_strand_mutations_C>G_c3'],
            'G-strand': ['g_strand_mutations_G>C_g1', 'g_strand_mutations_G>C_g2', 'g_strand_mutations_G>C_g3']
        },
        'C>T': {
            'C-strand': ['c_strand_mutations_C>T_c1', 'c_strand_mutations_C>T_c2', 'c_strand_mutations_C>T_c3'],
            'G-strand': ['g_strand_mutations_G>A_g1', 'g_strand_mutations_G>A_g2', 'g_strand_mutations_G>A_g3']
        },
        'G>A': {
            'G-strand': ['g_strand_mutations_G>A_g1', 'g_strand_mutations_G>A_g2', 'g_strand_mutations_G>A_g3'],
            'C-strand': ['c_strand_mutations_C>G_c1', 'c_strand_mutations_C>G_c2', 'c_strand_mutations_C>G_c3']
        },
        'G>C': {
            'G-strand': ['g_strand_mutations_G>C_g1', 'g_strand_mutations_G>C_g2', 'g_strand_mutations_G>C_g3'],
            'C-strand': ['c_strand_mutations_C>T_c1', 'c_strand_mutations_C>T_c2', 'c_strand_mutations_C>T_c3']
        },
        'G>T': {
            'G-strand': ['g_strand_mutations_G>T_g1', 'g_strand_mutations_G>T_g2', 'g_strand_mutations_G>T_g3'],
            'C-strand': ['c_strand_mutations_C>A_c1', 'c_strand_mutations_C>A_c2', 'c_strand_mutations_C>A_c3']
        }
    }
    # Ensure plots directory exists
    os.makedirs('plots', exist_ok=True)
    for idx, row in df.iterrows():
        filename = str(row['FileName']) if 'FileName' in row else f'sample_{idx}'
        # Remove file extension if present
        filename_base = os.path.splitext(filename)[0]
        output_path = os.path.join('plots', f'{filename_base}.png')
        plot_mutational_signature_row(row, mutation_types, mutation_columns, output_path)

def plot_spearman_with_age(csv_path):
    """
    For each numeric column in the CSV (except 'Age'), plot a scatter plot with Age on the x-axis and the column on the y-axis,
    calculate and display the Spearman correlation, and save each plot in 'spearman's plots' directory.
    Also outputs a CSV table with the r and p values for each column.
    """
    import scipy.stats as stats
    # Set seaborn style for consistency
    sns.set_theme(style="whitegrid", font_scale=1.1)
    # Read data
    df = pd.read_csv(csv_path)
    # Ensure output directory exists
    output_dir = "spearman's plots"
    os.makedirs(output_dir, exist_ok=True)
    # Only keep columns with numeric data and drop rows with missing Age
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if 'Age' not in numeric_cols:
        print("No 'Age' column found in the data.")
        return
    numeric_cols = [col for col in numeric_cols if col != 'Age']
    # Drop rows with missing Age
    df = df.dropna(subset=['Age'])
    spearman_results = []
    for col in numeric_cols:
        # Drop rows with missing values in the current column
        sub_df = df.dropna(subset=[col])
        if sub_df.shape[0] < 2:
            continue  # Not enough data to plot
        x = sub_df['Age']
        y = sub_df[col]
        # Calculate Spearman correlation
        corr, pval = stats.spearmanr(x, y)
        spearman_results.append({'Column': col, 'Spearman_r': corr, 'p_value': pval})
        # Plot using seaborn
        plt.figure(figsize=(8, 6))
        ax = sns.scatterplot(x=x, y=y)
        # Add seaborn regression (trendline) with no confidence interval
        if len(x) > 1:
            sns.regplot(x=x, y=y, scatter=False, ci=None, line_kws={'color': 'red', 'linestyle': '--'}, ax=ax)
        ax.set_xlabel('Age (years)', fontsize=12)
        ax.set_ylabel(col, fontsize=12)
        ax.set_title(f"Spearman's ρ = {corr:.2f} (p={pval:.2g})\n{col} vs Age", fontsize=14)
        plt.tight_layout()
        # Save plot
        safe_col = col.replace('/', '_').replace(' ', '_').replace('>', 'to').replace('<', 'lt').replace(':', '_')
        output_path = os.path.join(output_dir, f"{safe_col}_vs_Age.png")
        plt.savefig(output_path, dpi=200)
        plt.close()
    # Output CSV table of results
    results_df = pd.DataFrame(spearman_results)
    results_csv_path = os.path.join(output_dir, "spearman_results.csv")
    results_df.to_csv(results_csv_path, index=False)

# --- Composite Score Plotting ---
def plot_composite_score(csv_path, target_col='Age'):
    """
    Calculate a composite score from selected columns, plot it against the target column,
    and display the Spearman correlation.
    """
    import scipy.stats as stats
    # Columns to combine
    top_features = [
        'total_mutations_over_total_g_strand_2xrepeats_per_1k',
        'g_strand_mutations_A>C_a1_per_1k',
        'g_strand_mutations_T>G_t2_per_1k',
        'g_strand_mutations_T>G_t1_per_1k',
        'g_strand_mutations_G>A_g3_per_1k'
    ]
    df = pd.read_csv(csv_path)
    # Drop rows with missing target
    df = df.dropna(subset=[target_col])
    # Calculate composite score
    df['composite_score'] = df[top_features].mean(axis=1)
    # Calculate Spearman correlation
    r, p = stats.spearmanr(df['composite_score'], df[target_col])
    # Plot
    plt.figure(figsize=(8, 6))
    ax = sns.scatterplot(x=df[target_col], y=df['composite_score'])
    sns.regplot(x=df[target_col], y=df['composite_score'], scatter=False, ci=None, line_kws={'color': 'red', 'linestyle': '--'}, ax=ax)
    ax.set_xlabel(target_col, fontsize=12)
    ax.set_ylabel('Composite Score', fontsize=12)
    ax.set_title(f"Composite Score vs {target_col}\nSpearman's ρ = {r:.2f} (p={p:.2g})", fontsize=14)
    plt.tight_layout()
    # Save plot
    output_dir = "spearman's plots"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"composite_score_vs_{target_col}.png")
    plt.savefig(output_path, dpi=200)
    plt.close()
    print(f"Composite score plot saved as {output_path}")


def plot_composite_score_main():
    plot_composite_score('telomere_analysis.csv', target_col='Age')

def plot_mutational_signatures_main():
    plot_mutational_signatures('telomere_analysis.csv')
    print("Mutational signature plots saved in 'plots/' directory")

def plot_spearman_with_age_main():
    plot_spearman_with_age('telomere_analysis.csv')
    print("Spearman plots saved in 'spearman's plots/' directory")

if __name__ == "__main__":
    plot_mutational_signatures_main()
    plot_spearman_with_age_main()
    plot_composite_score_main()
