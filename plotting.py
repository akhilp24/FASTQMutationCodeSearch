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

    # For demonstration, we will use the G>T and G>A as C>T and C>A, etc.
    # You may need to adjust this mapping based on your actual data columns
    # Here, we assume G>T_g1/g2/g3 are C>A, G>A_g1/g2/g3 are C>T, etc.
    # This is a placeholder and should be mapped to your actual mutation context
    mutation_columns = [
        ['G_T_g1', 'G_T_g2', 'G_T_g3'],  # C>A
        [],                              # C>G (not present in your data)
        ['G_A_g1', 'G_A_g2', 'G_A_g3'],  # C>T
        [],                              # T>A (not present in your data)
        [],                              # T>C (not present in your data)
        [],                              # T>G (not present in your data)
    ]

    # Aggregate counts for each mutation type
    bar_heights = []
    bar_colors = []
    bar_labels = []
    total_mutations = df[sum(mutation_columns, [])].sum().sum()  # Total count of all mutations

    for i, (mut_label, color) in enumerate(mutation_types):
        cols = mutation_columns[i]
        if cols:
            values = df[cols].sum(axis=0)
            percentages = (values / total_mutations) * 100  # Convert counts to percentages
            bar_heights.extend(percentages)
            bar_colors.extend([color] * len(values))
            bar_labels.extend([mut_label] * len(values))
        else:
            # Add placeholders for missing types
            bar_heights.extend([0, 0, 0])
            bar_colors.extend([color] * 3)
            bar_labels.extend([mut_label] * 3)

    x = np.arange(len(bar_heights))

    fig, ax = plt.subplots(figsize=(6, 6))
    bars = ax.bar(x, bar_heights, color=bar_colors, edgecolor='black')

    # Add colored boxes and mutation type labels below each group
    group_positions = np.arange(0, len(bar_heights), 3)
    for i, (mut_label, color) in enumerate(mutation_types):
        # Draw colored box below bars
        # ax.add_patch(plt.Rectangle((group_positions[i] - 0.5, -max(bar_heights)*0.08), 3, max(bar_heights)*0.07, color=color, clip_on=False, zorder=2))
        # Draw label below colored box
        ax.text(group_positions[i] + 1, -max(bar_heights)*0.15, mut_label, ha='center', va='center', color='black' if color != 'black' else 'black', fontsize=16, fontweight='normal')

    ax.set_xticks([])
    ax.set_yticks(np.linspace(0, max(bar_heights), 5))
    ax.set_xlim(-0.5, len(bar_heights) - 0.5)
    ax.set_ylim(-max(bar_heights)*0.2, max(bar_heights) + max(bar_heights)*0.2)
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
