import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_mutational_signatures(csv_path, output_path='mutational_signatures.png'):
    """
    Plot mutational signatures (C>A, C>G, C>T, T>A, T>C, T>G) from a summary CSV file.
    Args:
        csv_path (str): Path to the summary CSV file.
        output_path (str): Path to save the generated plot.
    """
    df = pd.read_csv(csv_path)

    # Define the 6 mutation types and their colors
    mutation_types = [
        ('C>A', '#3DB7E9'),
        ('C>G', '#000000'),
        ('C>T', '#E71D25'),
        ('T>A', '#BEBEBE'),
        ('T>C', '#94C849'),
        ('T>G', '#F6C4C4'),
    ]

    # Aggregate columns for each mutation type (sum over all positions and both G- and C-strand contexts)
    mutation_columns = {
        'C>A': ['C_A_c1', 'C_A_c2', 'C_A_c3'],
        'C>G': ['C_G_c1', 'C_G_c2', 'C_G_c3'],
        'C>T': ['C_T_c1', 'C_T_c2', 'C_T_c3'],
        'T>A': ['T_A_t1', 'T_A_t2', 'T_A_t3'],
        'T>C': ['T_C_t1', 'T_C_t2', 'T_C_t3'],
        'T>G': ['T_G_t1', 'T_G_t2', 'T_G_t3'],
    }
    # Add G-strand context for C>X mutations (G>X_g1, g2, g3)
    mutation_columns['C>A'] += ['G_T_g1', 'G_T_g2', 'G_T_g3']  # G>T is C>A on the opposite strand
    mutation_columns['C>G'] += ['G_C_g1', 'G_C_g2', 'G_C_g3']  # G>C is C>G on the opposite strand
    mutation_columns['C>T'] += ['G_A_g1', 'G_A_g2', 'G_A_g3']  # G>A is C>T on the opposite strand

    # Aggregate counts for each mutation type
    bar_heights = []
    bar_colors = []
    bar_labels = []
    total_mutations = df[sum(mutation_columns.values(), [])].sum().sum()  # Total count of all mutations

    for mut_label, color in mutation_types:
        cols = mutation_columns[mut_label]
        values = df[cols].sum(axis=0)
        percentage = (values.sum() / total_mutations) * 100 if total_mutations > 0 else 0
        bar_heights.append(percentage)
        bar_colors.append(color)
        bar_labels.append(mut_label)

    x = np.arange(len(bar_heights))

    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.bar(x, bar_heights, color=bar_colors, edgecolor='black')

    # Add mutation type labels below each bar
    for i, (mut_label, color) in enumerate(mutation_types):
        ax.text(i, -max(bar_heights)*0.05, mut_label, ha='center', va='center', color='black', fontsize=14, fontweight='normal')

    ax.set_xticks([])
    ax.set_yticks(np.linspace(0, max(bar_heights), 5))
    ax.set_xlim(-0.5, len(bar_heights) - 0.5)
    ax.set_ylim(-max(bar_heights)*0.1, max(bar_heights) + max(bar_heights)*0.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_ylabel('Percentage of Single Base Modifications', fontsize=14)
    ax.set_title('Mutational Signatures', fontsize=18, fontweight='bold', pad=30)
    ax.yaxis.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python plotting.py <csv_path> [output_path]")
        sys.exit(1)
    csv_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else 'mutational_signatures.png'
    plot_mutational_signatures(csv_path, output_path)
