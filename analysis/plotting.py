import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns

import numbers

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
        'total_mutations_over_total_g_strand_3xrepeats_per_1k',
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

def plot_mutation_r_heatmap(csv_path, target_col='Age'):
    """
    Plot a heatmap of Spearman r values between each normalized mutation column (per_1k/per1k) and the target column (e.g., Age),
    as well as the total mutation count column if present. Also plot a clustered heatmap of these values for all samples.
    """
    import scipy.stats as stats
    # Read data
    df = pd.read_csv(csv_path)
    # Only keep columns with numeric data and drop rows with missing target
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if target_col not in numeric_cols:
        print(f"No '{target_col}' column found in the data.")
        return
    # Only use columns with 'per_1k' or 'per1k' in the name, plus total mutation count if present
    mutation_cols = [col for col in numeric_cols if ('per_1k' in col or 'per1k' in col)]
    total_mut_col = None
    for col in numeric_cols:
        if 'total_mutation' in col and col not in mutation_cols:
            total_mut_col = col
            break
    if total_mut_col:
        mutation_cols.append(total_mut_col)
    if not mutation_cols:
        print("No normalized 'per_1k' or 'per1k' mutation columns found.")
        return
    df = df.dropna(subset=[target_col])
    # Calculate Spearman r for each mutation column
    r_values = []
    for col in mutation_cols:
        sub_df = df.dropna(subset=[col])
        if sub_df.shape[0] < 2:
            r = float('nan')
        else:
            r, _ = stats.spearmanr(sub_df[col], sub_df[target_col])
        r_values.append(r)
    # Create DataFrame for heatmap
    r_df = pd.DataFrame({'Mutation': mutation_cols, 'Spearman_r': r_values})
    r_df = r_df.set_index('Mutation')
    # Plot heatmap of r values
    plt.figure(figsize=(max(8, len(mutation_cols) * 0.4), 2.5))
    sns.heatmap(r_df.T, annot=True, cmap='coolwarm', center=0, cbar_kws={'label': "Spearman's r"})
    plt.title(f"Spearman r values: Normalized Mutations vs {target_col}")
    plt.yticks(rotation=0)
    plt.tight_layout()
    output_dir = "spearman's plots"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"mutation_r_heatmap_vs_{target_col}.png")
    plt.savefig(output_path, dpi=200)
    plt.close()
    print(f"Mutation r heatmap saved as {output_path}")

    # Clustered heatmap of normalized mutation values for all samples
    # (Removed as per user request)
    # if len(mutation_cols) > 1:
    #     import warnings
    #     with warnings.catch_warnings():
    #         warnings.simplefilter("ignore")
    #         sns.clustermap(df[mutation_cols], cmap='viridis', standard_scale=1, yticklabels=False)
    #     plt.suptitle('Clustered Heatmap of Normalized Mutation Values (Z-score by column)', y=1.02)
    #     plt.tight_layout()
    #     output_path2 = os.path.join(output_dir, f"mutation_value_clustermap.png")
    #     plt.savefig(output_path2, dpi=200)
    #     plt.close()
    #     print(f"Clustered mutation value heatmap saved as {output_path2}")


def plot_mutation_r_heatmap_main():
    plot_mutation_r_heatmap('telomere_analysis.csv', target_col='Age')

def plot_pairwise_r_heatmap(csv_path):
    """
    Plot a heatmap of pairwise Spearman r values between all relevant columns (per_1k/per1k, total mutation count, Age, and telomere columns).
    The diagonal is masked (crossed out).
    """
    import scipy.stats as stats
    # Read data
    df = pd.read_csv(csv_path)
    # Only use columns with 'per_1k' or 'per1k' in the name, plus total mutation count if present
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    mutation_cols = [col for col in numeric_cols if ('per_1k' in col or 'per1k' in col)]
    total_mut_col = None
    for col in numeric_cols:
        if 'total_mutation' in col and col not in mutation_cols:
            total_mut_col = col
            break
    if total_mut_col:
        mutation_cols.append(total_mut_col)
    # Add Age if present
    if 'Age' in df.columns and 'Age' in numeric_cols and 'Age' not in mutation_cols:
        mutation_cols.append('Age')
    # Add telomere columns (case-insensitive)
    telomere_cols = [col for col in df.columns if 'telomere' in col.lower() and col in numeric_cols and col not in mutation_cols]
    mutation_cols.extend(telomere_cols)
    if not mutation_cols:
        print("No relevant columns found for pairwise heatmap.")
        return
    # Compute pairwise Spearman r matrix
    n = len(mutation_cols)
    r_matrix = np.zeros((n, n))
    for i, col1 in enumerate(mutation_cols):
        for j, col2 in enumerate(mutation_cols):
            # Drop rows with missing values in either column
            sub_df = df[[col1, col2]].dropna()
            if sub_df.shape[0] < 2:
                r = np.nan
            else:
                r_val = stats.spearmanr(sub_df[col1], sub_df[col2])[0]
                # Ensure r_val is a scalar float
                if isinstance(r_val, numbers.Number) and not isinstance(r_val, (list, tuple, np.ndarray)):
                    r = float(r_val)
                else:
                    r = np.nan
            r_matrix[i, j] = r
    # Mask the diagonal
    mask = np.eye(n, dtype=bool)
    # Plot heatmap with improved readability
    fig_width = max(10, n * 0.7)
    fig_height = max(8, n * 0.7)
    plt.figure(figsize=(fig_width, fig_height))
    ax = sns.heatmap(
        r_matrix,
        annot=True,
        fmt=".2f",
        cmap="coolwarm",
        center=0,
        mask=mask,
        xticklabels=mutation_cols,
        yticklabels=mutation_cols,
        cbar_kws={'label': "Spearman's r"},
        annot_kws={"size": 8}
    )
    # Cross out the diagonal
    for i in range(n):
        ax.add_patch(plt.Rectangle((i, i), 1, 1, fill=False, edgecolor='black', lw=2, hatch='xx'))
    plt.title("Pairwise Spearman r Heatmap (Normalized Mutations, Age, Telomere)", fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yticks(fontsize=9)
    plt.tight_layout()
    output_dir = "spearman's plots"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "pairwise_mutation_r_heatmap.png")
    plt.savefig(output_path, dpi=200)
    plt.close()
    print(f"Pairwise mutation r heatmap saved as {output_path}")

def plot_pairwise_r_heatmap_main():
    plot_pairwise_r_heatmap('telomere_analysis.csv')

if __name__ == "__main__":
    plot_mutational_signatures_main()
    plot_spearman_with_age_main()
    plot_composite_score_main()
    plot_mutation_r_heatmap_main()
    plot_pairwise_r_heatmap_main()
