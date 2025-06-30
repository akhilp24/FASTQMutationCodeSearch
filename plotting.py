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
            'C-strand': ['C_A_c1', 'C_A_c2', 'C_A_c3'],
            'G-strand': ['G_T_g1', 'G_T_g2', 'G_T_g3']
        },
        'C>G': {
            'C-strand': ['C_G_c1', 'C_G_c2', 'C_G_c3'],
            'G-strand': ['G_C_g1', 'G_C_g2', 'G_C_g3']
        },
        'C>T': {
            'C-strand': ['C_T_c1', 'C_T_c2', 'C_T_c3'],
            'G-strand': ['G_A_g1', 'G_A_g2', 'G_A_g3']
        },
        'G>A': {
            'G-strand': ['G_A_g1', 'G_A_g2', 'G_A_g3'],
            'C-strand': ['C_G_c1', 'C_G_c2', 'C_G_c3']
        },
        'G>C': {
            'G-strand': ['G_C_g1', 'G_C_g2', 'G_C_g3'],
            'C-strand': ['C_T_c1', 'C_T_c2', 'C_T_c3']
        },
        'G>T': {
            'G-strand': ['G_T_g1', 'G_T_g2', 'G_T_g3'],
            'C-strand': ['C_A_c1', 'C_A_c2', 'C_A_c3']
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

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python plotting.py <csv_path>")
        sys.exit(1)
    csv_path = sys.argv[1]
    plot_mutational_signatures(csv_path)
