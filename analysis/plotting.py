import json
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns
from scipy.optimize import curve_fit
from scipy import stats

import numbers


def _get_patterns_version(patterns_file_path):
    """Read version string from patterns JSON (path passed from main.py)."""
    try:
        with open(patterns_file_path) as f:
            return json.load(f).get('version', 'unknown')
    except Exception:
        return 'unknown'


def plot_mutational_signature_row(row, mutation_types, mutation_columns, output_path, version):
    # Set seaborn style for better-looking plots
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    
    # Aggregate counts for each mutation type, position, and strand context for a single row
    bar_heights = []
    bar_colors = []
    bar_labels = []
    # Get all columns for total calculation - only use columns that exist in the data
    all_columns = []
    for mut_type, contexts in mutation_columns.items():
        for context, cols in contexts.items():
            # Filter to only include columns that exist in the row
            existing_cols = [col for col in cols if col in row.index]
            all_columns.extend(existing_cols)
    
    if not all_columns:
        print(f"Warning: No valid mutation columns found for {row.get('FileName', 'unknown')}")
        return
    
    total_mutations = row[all_columns].sum()
    
    for mut_label, color in mutation_types:
        contexts = mutation_columns[mut_label]
        for context_name, cols in contexts.items():
            for i, col in enumerate(cols):
                if col in row.index:  # Only process existing columns
                    value = row[col]
                    percentage = (value / total_mutations) * 100 if total_mutations > 0 else 0
                    bar_heights.append(percentage)
                    bar_colors.append(color)
                    bar_labels.append(f"{mut_label} {context_name} pos{i+1}")
    
    if not bar_heights:
        print(f"Warning: No valid data found for {row.get('FileName', 'unknown')}")
        return
    
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
    title = f"Mutational Signatures by Position and Strand Context [{version}]\nFile: {filename} | Age: {age} years"
    ax.set_title(title, fontsize=18, fontweight='bold', pad=30)
    
    # Improve layout
    plt.tight_layout(rect=[0, 0.15, 1, 0.95])
    
    # Save with high quality
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def plot_mutational_signatures(csv_path, patterns_file_path):
    # Set global seaborn style
    sns.set_theme(style="whitegrid", font_scale=1.1)
    version = _get_patterns_version(patterns_file_path)
    
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
        },
       
    }
    # Ensure plots directory exists
    os.makedirs('plots', exist_ok=True)
    for idx, row in df.iterrows():
        filename = str(row['FileName']) if 'FileName' in row else f'sample_{idx}'
        filename_base = os.path.splitext(filename)[0]
        output_path = os.path.join('plots', f'{filename_base}.png')
        plot_mutational_signature_row(row, mutation_types, mutation_columns, output_path, version)

def plot_spearman_with_age(csv_path, patterns_file_path):
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
        version = _get_patterns_version(patterns_file_path)
        ax.set_title(f"Spearman's ρ = {corr:.2f} (p={pval:.2g})\n{col} vs Age [{version}]", fontsize=14)
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
def plot_composite_score(csv_path, target_col='Age', patterns_file_path=None):
    """
    Calculate a composite score from selected columns, plot it against the target column,
    and display the Spearman correlation.
    """
    import scipy.stats as stats
    # Columns to combine - updated to use read-based metrics (frameshifts disabled)
    top_features = [
        'total_mutations_per_1k_reads',  # Updated to use read-based normalization
        'g_strand_mutations_A>C_a1_per_1k',
        'g_strand_mutations_T>G_t2_per_1k',
        'g_strand_mutations_T>G_t1_per_1k',
        'g_strand_mutations_G>A_g3_per_1k',
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
    version = _get_patterns_version(patterns_file_path) if patterns_file_path else 'unknown'
    ax.set_title(f"Composite Score vs {target_col} [{version}]\nSpearman's ρ = {r:.2f} (p={p:.2g})", fontsize=14)
    plt.tight_layout()
    # Save plot
    output_dir = "spearman's plots"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"composite_score_vs_{target_col}.png")
    plt.savefig(output_path, dpi=200)
    plt.close()
    print(f"Composite score plot saved as {output_path}")


def plot_composite_score_main(patterns_file_path):
    plot_composite_score('telomere_analysis.csv', target_col='Age', patterns_file_path=patterns_file_path)

def plot_mutational_signatures_main(patterns_file_path):
    plot_mutational_signatures('telomere_analysis.csv', patterns_file_path)
    print("Mutational signature plots saved in 'plots/' directory")

def plot_spearman_with_age_main(patterns_file_path):
    plot_spearman_with_age('telomere_analysis.csv', patterns_file_path)
    print("Spearman plots saved in 'spearman's plots/' directory")

def plot_mutation_r_heatmap(csv_path, target_col='Age', patterns_file_path=None):
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
    
    # Add the new read-based total mutations column
    if 'total_mutations_per_1k_reads' in df.columns and 'total_mutations_per_1k_reads' not in mutation_cols:
        mutation_cols.append('total_mutations_per_1k_reads')
    
    # Add legacy total mutation columns for backward compatibility
    total_mut_col = None
    for col in numeric_cols:
        if 'total_mutation' in col and col not in mutation_cols:
            total_mut_col = col
            break
    if total_mut_col:
        mutation_cols.append(total_mut_col)
        
    # Add frameshift-specific columns
    frameshift_cols = [col for col in numeric_cols if 'frameshift' in col.lower() and 'per_1k' in col]
    mutation_cols.extend([col for col in frameshift_cols if col not in mutation_cols])
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
    version = _get_patterns_version(patterns_file_path) if patterns_file_path else 'unknown'
    plt.title(f"Spearman r values: Normalized Mutations vs {target_col} [{version}]")
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


def plot_mutation_r_heatmap_main(patterns_file_path):
    plot_mutation_r_heatmap('telomere_analysis.csv', target_col='Age', patterns_file_path=patterns_file_path)

def plot_pairwise_r_heatmap(csv_path, patterns_file_path=None):
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
    
    # Add the new read-based total mutations column
    if 'total_mutations_per_1k_reads' in df.columns and 'total_mutations_per_1k_reads' not in mutation_cols:
        mutation_cols.append('total_mutations_per_1k_reads')
    
    # Add legacy total mutation columns for backward compatibility
    total_mut_col = None
    for col in numeric_cols:
        if 'total_mutation' in col and col not in mutation_cols:
            total_mut_col = col
            break
    if total_mut_col:
        mutation_cols.append(total_mut_col)
        
    # Add frameshift-specific columns
    frameshift_cols = [col for col in numeric_cols if 'frameshift' in col.lower() and 'per_1k' in col]
    mutation_cols.extend([col for col in frameshift_cols if col not in mutation_cols])
    
    # Add Total_Reads column if present
    if 'Total_Reads' in df.columns and 'Total_Reads' in numeric_cols and 'Total_Reads' not in mutation_cols:
        mutation_cols.append('Total_Reads')
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
    version = _get_patterns_version(patterns_file_path) if patterns_file_path else 'unknown'
    plt.title(f"Pairwise Spearman r Heatmap (Normalized Mutations, Age, Telomere) [{version}]", fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yticks(fontsize=9)
    plt.tight_layout()
    output_dir = "spearman's plots"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "pairwise_mutation_r_heatmap.png")
    plt.savefig(output_path, dpi=200)
    plt.close()
    print(f"Pairwise mutation r heatmap saved as {output_path}")

def plot_pairwise_r_heatmap_main(patterns_file_path):
    plot_pairwise_r_heatmap('telomere_analysis.csv', patterns_file_path=patterns_file_path)

def curve_fitting_analysis(csv_path, output_dir="curve_fitting_plots", patterns_file_path=None):
    """
    Perform curve fitting analysis for telomere length vs age and mutation rate vs age.
    Try multiple curve types (linear, exponential, logarithmic, polynomial, power) and export the best-fit plots.
    
    The function fits 5 different curve types to each variable and determines the best fit based on R-squared values:
    - Linear: y = ax + b
    - Exponential: y = ae^(bx) + c  
    - Logarithmic: y = a*ln(x+1) + b
    - Polynomial (cubic): y = ax³ + bx² + cx + d
    - Power: y = a(x+1)^b + c
    
    For each variable, generates a plot showing all fitted curves with the best fit highlighted,
    identifies outlier points (>2 standard deviations from best-fit curve), and annotates them
    with sample names. Saves detailed results including R-squared values and fitted parameters to CSV.
    
    Args:
        csv_path: Path to CSV file with Age, Telomere_Length, and mutation rate columns
        output_dir: Directory to save output plots and results (default: "curve_fitting_plots")
    
    Outputs:
        - Individual curve fitting plots for each variable with outliers highlighted and annotated
        - curve_fitting_results.csv with R-squared values and parameters for all fits
        - Console summary of best-fit curves and outlier samples for each variable
    """
    # Define curve fitting functions
    def linear_func(x, a, b):
        return a * x + b
    
    def exponential_func(x, a, b, c):
        return a * np.exp(b * x) + c
    
    def logarithmic_func(x, a, b):
        return a * np.log(x + 1) + b  # +1 to avoid log(0)
    
    def polynomial_func(x, a, b, c, d):
        return a * x**3 + b * x**2 + c * x + d
    
    def power_func(x, a, b, c):
        return a * (x + 1)**b + c  # +1 to avoid negative bases
    
    # Read data
    df = pd.read_csv(csv_path)
    
    # Drop rows with missing Age
    df = df.dropna(subset=['Age'])
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Set seaborn style
    sns.set_theme(style="whitegrid", font_scale=1.1)
    
    results = []
    
    # 1. Telomere Length vs Age Analysis
    if 'Telomere_Length' in df.columns:
        telomere_df = df.dropna(subset=['Telomere_Length'])
        if len(telomere_df) >= 4:  # Need at least 4 points for fitting
            x_data = telomere_df['Age'].values
            y_data = telomere_df['Telomere_Length'].values
            
            # Try different curve types
            curve_types = [
                ('Linear', linear_func, [1, 1]),
                ('Exponential', exponential_func, [1, 0.01, 1]),
                ('Logarithmic', logarithmic_func, [1, 1]),
                ('Polynomial', polynomial_func, [0.01, 0.1, 1, 1]),
                ('Power', power_func, [1, 0.5, 1])
            ]
            
            best_fit = None
            best_r_squared = -np.inf
            
            plt.figure(figsize=(14, 10))  # Larger figure to accommodate annotations
            colors = ['red', 'green', 'blue', 'orange', 'purple']
            
            for i, (name, func, initial_guess) in enumerate(curve_types):
                try:
                    # Perform curve fitting
                    popt, pcov = curve_fit(func, x_data, y_data, p0=initial_guess, maxfev=5000)
                    
                    # Calculate R-squared
                    y_pred = func(x_data, *popt)
                    ss_res = np.sum((y_data - y_pred) ** 2)
                    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
                    r_squared = 1 - (ss_res / ss_tot)
                    
                    # Store results
                    results.append({
                        'Variable': 'Telomere_Length',
                        'Curve_Type': name,
                        'R_squared': r_squared,
                        'Parameters': popt.tolist(),
                        'Parameter_Errors': np.sqrt(np.diag(pcov)).tolist()
                    })
                    
                    # Plot the fitted curve
                    x_smooth = np.linspace(x_data.min(), x_data.max(), 100)
                    y_smooth = func(x_smooth, *popt)
                    plt.plot(x_smooth, y_smooth, color=colors[i], linestyle='--', 
                            label=f'{name} (R² = {r_squared:.3f})', linewidth=2)
                    
                    # Track best fit
                    if r_squared > best_r_squared:
                        best_r_squared = r_squared
                        best_fit = (name, func, popt)
                        
                except Exception as e:
                    print(f"Could not fit {name} curve to Telomere_Length vs Age: {e}")
            
            # Plot original data
            plt.scatter(x_data, y_data, alpha=0.7, s=50, color='black', label='Data points')
            
            # Highlight best fit and identify outliers
            if best_fit:
                name, func, popt = best_fit
                x_smooth = np.linspace(x_data.min(), x_data.max(), 100)
                y_smooth = func(x_smooth, *popt)
                plt.plot(x_smooth, y_smooth, color='red', linewidth=3, 
                        label=f'Best fit: {name} (R² = {best_r_squared:.3f})')
                
                # Identify and annotate outliers
                y_pred = func(x_data, *popt)
                residuals = y_data - y_pred
                residual_std = np.std(residuals)
                outlier_threshold = 2 * residual_std  # Points > 2 std deviations
                
                outlier_indices = np.where(np.abs(residuals) > outlier_threshold)[0]
                if len(outlier_indices) > 0:
                    # Plot outliers with different color
                    plt.scatter(x_data[outlier_indices], y_data[outlier_indices], 
                              alpha=0.9, s=80, color='orange', edgecolor='red', linewidth=2,
                              label=f'Outliers (>{outlier_threshold:.1f} from curve)', zorder=5)
                    
                    print(f"\n--- Telomere Length Outliers ---")
                    print(f"Outlier threshold: ±{outlier_threshold:.1f} bp")
                    
                    # Annotate outliers with sample names
                    for idx in outlier_indices:
                        sample_name = telomere_df.iloc[idx]['FileName'] if 'FileName' in telomere_df.columns else f'Sample_{idx}'
                        age = x_data[idx]
                        telomere_length = y_data[idx]
                        residual = residuals[idx]
                        print(f"  {sample_name}: Age={age:.1f}, TL={telomere_length:.1f}bp, Residual={residual:+.1f}bp")
                        
                        # Shorten sample name for plot if too long
                        display_name = sample_name
                        if len(str(sample_name)) > 15:
                            display_name = str(sample_name)[:12] + '...'
                        plt.annotate(display_name, 
                                   (x_data[idx], y_data[idx]),
                                   xytext=(10, 10), textcoords='offset points',
                                   fontsize=8, ha='left', va='bottom',
                                   bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                                   arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.1'))
            
            plt.xlabel('Age (years)', fontsize=12)
            plt.ylabel('Telomere Length (bp)', fontsize=12)
            version = _get_patterns_version(patterns_file_path) if patterns_file_path else 'unknown'
            plt.title(f'Curve Fitting: Telomere Length vs Age [{version}]', fontsize=14, fontweight='bold')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            
            # Save plot
            output_path = os.path.join(output_dir, 'telomere_length_vs_age_curve_fitting.png')
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Telomere length curve fitting plot saved: {output_path}")
    
    # 2. Mutation Rate vs Age Analysis
    # Find mutation rate columns (normalized per 1k)
    mutation_rate_cols = [col for col in df.columns if 'per_1k' in col and col != 'Telomere_Length']
    
    if mutation_rate_cols:
        # Use total mutation rate or composite score if available - prioritize read-based metrics
        target_cols = []
        for col in ['total_mutations_per_1k_reads', 'mutation_rate_normalized_by_length', 'composite_score', 
                   # 'total_g_strand_frameshifts_per_1k', 'total_c_strand_frameshifts_per_1k',
                   'g_strand_mutations_sum_per_1k', 'c_strand_mutations_sum_per_1k']:
            if col in df.columns:
                target_cols.append(col)
        
        # If no composite columns, use first few individual mutation rate columns
        if not target_cols:
            target_cols = mutation_rate_cols[:3]  # Take first 3 columns
        
        for col in target_cols:
            mutation_df = df.dropna(subset=[col])
            if len(mutation_df) >= 4:  # Need at least 4 points for fitting
                x_data = mutation_df['Age'].values
                y_data = mutation_df[col].values
                
                # Try different curve types
                curve_types = [
                    ('Linear', linear_func, [1, 1]),
                    ('Exponential', exponential_func, [1, 0.01, 1]),
                    ('Logarithmic', logarithmic_func, [1, 1]),
                    ('Polynomial', polynomial_func, [0.01, 0.1, 1, 1]),
                    ('Power', power_func, [1, 0.5, 1])
                ]
                
                best_fit = None
                best_r_squared = -np.inf
                
                plt.figure(figsize=(14, 10))  # Larger figure to accommodate annotations
                colors = ['red', 'green', 'blue', 'orange', 'purple']
                
                for i, (name, func, initial_guess) in enumerate(curve_types):
                    try:
                        # Perform curve fitting
                        popt, pcov = curve_fit(func, x_data, y_data, p0=initial_guess, maxfev=5000)
                        
                        # Calculate R-squared
                        y_pred = func(x_data, *popt)
                        ss_res = np.sum((y_data - y_pred) ** 2)
                        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
                        r_squared = 1 - (ss_res / ss_tot)
                        
                        # Store results
                        results.append({
                            'Variable': col,
                            'Curve_Type': name,
                            'R_squared': r_squared,
                            'Parameters': popt.tolist(),
                            'Parameter_Errors': np.sqrt(np.diag(pcov)).tolist()
                        })
                        
                        # Plot the fitted curve
                        x_smooth = np.linspace(x_data.min(), x_data.max(), 100)
                        y_smooth = func(x_smooth, *popt)
                        plt.plot(x_smooth, y_smooth, color=colors[i], linestyle='--', 
                                label=f'{name} (R² = {r_squared:.3f})', linewidth=2)
                        
                        # Track best fit
                        if r_squared > best_r_squared:
                            best_r_squared = r_squared
                            best_fit = (name, func, popt)
                            
                    except Exception as e:
                        print(f"Could not fit {name} curve to {col} vs Age: {e}")
                
                # Plot original data
                plt.scatter(x_data, y_data, alpha=0.7, s=50, color='black', label='Data points')
                
                # Highlight best fit and identify outliers
                if best_fit:
                    name, func, popt = best_fit
                    x_smooth = np.linspace(x_data.min(), x_data.max(), 100)
                    y_smooth = func(x_smooth, *popt)
                    plt.plot(x_smooth, y_smooth, color='red', linewidth=3, 
                            label=f'Best fit: {name} (R² = {best_r_squared:.3f})')
                    
                    # Identify and annotate outliers
                    y_pred = func(x_data, *popt)
                    residuals = y_data - y_pred
                    residual_std = np.std(residuals)
                    outlier_threshold = 2 * residual_std  # Points > 2 std deviations
                    
                    outlier_indices = np.where(np.abs(residuals) > outlier_threshold)[0]
                    if len(outlier_indices) > 0:
                        # Plot outliers with different color
                        plt.scatter(x_data[outlier_indices], y_data[outlier_indices], 
                                  alpha=0.9, s=80, color='orange', edgecolor='red', linewidth=2,
                                  label=f'Outliers (>{outlier_threshold:.3f} from curve)', zorder=5)
                        
                        print(f"\n--- {col.replace('_', ' ').title()} Outliers ---")
                        print(f"Outlier threshold: ±{outlier_threshold:.3f}")
                        
                        # Annotate outliers with sample names
                        for idx in outlier_indices:
                            sample_name = mutation_df.iloc[idx]['FileName'] if 'FileName' in mutation_df.columns else f'Sample_{idx}'
                            age = x_data[idx]
                            value = y_data[idx]
                            residual = residuals[idx]
                            print(f"  {sample_name}: Age={age:.1f}, Value={value:.3f}, Residual={residual:+.3f}")
                            
                            # Shorten sample name for plot if too long
                            display_name = sample_name
                            if len(str(sample_name)) > 15:
                                display_name = str(sample_name)[:12] + '...'
                            plt.annotate(display_name, 
                                       (x_data[idx], y_data[idx]),
                                       xytext=(10, 10), textcoords='offset points',
                                       fontsize=8, ha='left', va='bottom',
                                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.1'))
                
                plt.xlabel('Age (years)', fontsize=12)
                plt.ylabel(col.replace('_', ' ').title(), fontsize=12)
                version = _get_patterns_version(patterns_file_path) if patterns_file_path else 'unknown'
                plt.title(f'Curve Fitting: {col.replace("_", " ").title()} vs Age [{version}]', fontsize=14, fontweight='bold')
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                
                # Save plot
                safe_col = col.replace('/', '_').replace(' ', '_').replace('>', 'to').replace('<', 'lt').replace(':', '_')
                output_path = os.path.join(output_dir, f'{safe_col}_vs_age_curve_fitting.png')
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"Mutation rate curve fitting plot saved: {output_path}")
    
    # Save results to CSV
    if results:
        results_df = pd.DataFrame(results)
        results_csv_path = os.path.join(output_dir, "curve_fitting_results.csv")
        results_df.to_csv(results_csv_path, index=False)
        print(f"Curve fitting results saved: {results_csv_path}")
        
        # Print summary of best fits
        print("\n=== CURVE FITTING SUMMARY ===")
        for variable in results_df['Variable'].unique():
            var_results = results_df[results_df['Variable'] == variable]
            best_result = var_results.loc[var_results['R_squared'].idxmax()]
            print(f"{variable}:")
            print(f"  Best fit: {best_result['Curve_Type']} (R² = {best_result['R_squared']:.4f})")
            print(f"  Parameters: {best_result['Parameters']}")
            print()

def curve_fitting_analysis_main(patterns_file_path):
    """Main function to run curve fitting analysis"""
    curve_fitting_analysis('telomere_analysis.csv', patterns_file_path=patterns_file_path)


if __name__ == "__main__":
    _dir = os.path.dirname(__file__)
    _path = os.path.join(_dir, 'telomere_patterns_2x.json')
    plot_mutational_signatures_main(_path)
    plot_spearman_with_age_main(_path)
    plot_composite_score_main(_path)
    plot_mutation_r_heatmap_main(_path)
    plot_pairwise_r_heatmap_main(_path)
    curve_fitting_analysis_main(_path)
    # plot_frameshift_analysis_main()  # Frameshift analysis disabled
