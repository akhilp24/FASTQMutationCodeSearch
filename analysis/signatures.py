import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns

def load_telomere_data(file_path):
    """Load telomere analysis data from CSV file."""
    df = pd.read_csv(file_path)
    return df

def calculate_sbs_signatures(df):
    """
    Calculate SBS signature scores based on COSMIC mutational signatures.
    
    SBS1: Clock-like signature characterized by high C>T mutations
    SBS4: Tobacco smoking signature characterized by high C>A mutations  
    SBS5: Clock-like signature characterized by high C>T and T>C mutations
    SBS18: Reactive oxygen species signature characterized by high C>A and C>T mutations
    """
    signatures = {}
    
    # SBS1: C>T mutations (clock-like signature)
    # Sum all C>T mutations from both strands (C>T on C-strand + G>A on G-strand)
    ct_mutations = 0
    if 'c_strand_C>T_per_1k' in df.columns:
        ct_mutations += df['c_strand_C>T_per_1k']
    if 'g_strand_G>A_per_1k' in df.columns:
        ct_mutations += df['g_strand_G>A_per_1k']
    signatures['SBS1_CT_score'] = ct_mutations
    
    # SBS4: C>A mutations (tobacco smoking signature)  
    # Sum all C>A mutations from both strands (C>A on C-strand + G>T on G-strand)
    ca_mutations = 0
    if 'c_strand_C>A_per_1k' in df.columns:
        ca_mutations += df['c_strand_C>A_per_1k']
    if 'g_strand_G>T_per_1k' in df.columns:
        ca_mutations += df['g_strand_G>T_per_1k']
    signatures['SBS4_CA_score'] = ca_mutations
    
    # SBS5: C>T and T>C mutations (clock-like signature)
    # Sum C>T + T>C from both strands
    ct_tc_mutations = ct_mutations  # C>T already calculated above
    if 'c_strand_T>C_per_1k' in df.columns:
        ct_tc_mutations += df['c_strand_T>C_per_1k']
    if 'g_strand_A>G_per_1k' in df.columns:
        ct_tc_mutations += df['g_strand_A>G_per_1k']
    signatures['SBS5_CT_TC_score'] = ct_tc_mutations
    
    # SBS18: C>A and C>T mutations (reactive oxygen species signature)
    # Sum C>A + C>T from both strands
    ca_ct_mutations = ca_mutations + ct_mutations  # Both already calculated above
    signatures['SBS18_CA_CT_score'] = ca_ct_mutations
    
    return signatures

def calculate_correlations_with_age(df, signatures):
    """Calculate Spearman correlations between signature scores and age."""
    correlations = {}
    
    # Only use samples with age data
    age_mask = df['Age'].notna()
    age_data = df[age_mask]['Age']
    
    for sig_name, sig_scores in signatures.items():
        sig_data = sig_scores[age_mask]
        
        # Calculate Spearman correlation
        correlation, p_value = spearmanr(age_data, sig_data)
        
        correlations[sig_name] = {
            'correlation': correlation,
            'p_value': p_value,
            'n_samples': len(age_data)
        }
    
    return correlations

def rank_signatures(correlations):
    """Rank signatures by absolute correlation strength with age."""
    # Create list of tuples for sorting
    sig_rankings = []
    for sig_name, stats in correlations.items():
        sig_rankings.append((
            sig_name,
            abs(stats['correlation']),
            stats['correlation'],
            stats['p_value'],
            stats['n_samples']
        ))
    
    # Sort by absolute correlation strength (descending)
    sig_rankings.sort(key=lambda x: x[1], reverse=True)
    
    return sig_rankings

def create_visualization(df, signatures, correlations):
    """Create visualizations of signature correlations with age."""
    # Setup the plotting area
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('COSMIC SBS Mutational Signatures vs Age', fontsize=16, fontweight='bold')
    
    # Only use samples with age data
    age_mask = df['Age'].notna()
    age_data = df[age_mask]['Age']
    
    # List of signatures to plot
    sig_names = ['SBS1_CT_score', 'SBS4_CA_score', 'SBS5_CT_TC_score', 'SBS18_CA_CT_score']
    sig_titles = [
        'SBS1: C>T (Clock-like)',
        'SBS4: C>A (Tobacco smoking)', 
        'SBS5: C>T + T>C (Clock-like)',
        'SBS18: C>A + C>T (Reactive oxygen species)'
    ]
    
    # Create scatter plots
    for i, (sig_name, title) in enumerate(zip(sig_names, sig_titles)):
        row = i // 2
        col = i % 2
        ax = axes[row, col]
        
        sig_data = signatures[sig_name][age_mask]
        
        # Create scatter plot
        ax.scatter(age_data, sig_data, alpha=0.6, s=50)
        
        # Add trend line
        z = np.polyfit(age_data, sig_data, 1)
        p = np.poly1d(z)
        ax.plot(age_data, p(age_data), "r--", alpha=0.8)
        
        # Add correlation info to plot
        corr = correlations[sig_name]['correlation']
        p_val = correlations[sig_name]['p_value']
        ax.text(0.05, 0.95, f'r = {corr:.3f}\np = {p_val:.3e}', 
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax.set_xlabel('Age (years)')
        ax.set_ylabel('Mutations per 1k bases')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

def print_detailed_results(correlations, sig_rankings):
    """Print detailed results of the signature analysis."""
    print("="*80)
    print("COSMIC SBS Mutational Signatures Analysis")
    print("="*80)
    print()
    
    print("SIGNATURE DEFINITIONS:")
    print("-" * 40)
    print("• SBS1: Clock-like signature characterized by C>T mutations (spontaneous deamination of 5-methylcytosine)")
    print("• SBS4: Tobacco smoking signature characterized by C>A mutations")  
    print("• SBS5: Clock-like signature characterized by C>T and T>C mutations")
    print("• SBS18: Reactive oxygen species signature characterized by C>A and C>T mutations")
    print()
    
    print("SIGNATURE RANKINGS BY CORRELATION WITH AGE:")
    print("-" * 50)
    print(f"{'Rank':<4} {'Signature':<20} {'|r|':<8} {'r':<8} {'p-value':<12} {'N':<6}")
    print("-" * 50)
    
    for rank, (sig_name, abs_corr, corr, p_val, n_samples) in enumerate(sig_rankings, 1):
        # Clean up signature name for display
        display_name = sig_name.replace('_score', '').replace('_', ' ')
        print(f"{rank:<4} {display_name:<20} {abs_corr:<8.3f} {corr:<8.3f} {p_val:<12.3e} {n_samples:<6}")
    
    print()
    print("DETAILED STATISTICS:")
    print("-" * 30)
    for sig_name, stats in correlations.items():
        display_name = sig_name.replace('_score', '').replace('_', ' ')
        print(f"{display_name}:")
        print(f"  Spearman correlation: {stats['correlation']:.4f}")
        print(f"  P-value: {stats['p_value']:.4e}")
        print(f"  Sample size: {stats['n_samples']}")
        
        # Interpretation
        abs_corr = abs(stats['correlation'])
        if abs_corr >= 0.7:
            strength = "strong"
        elif abs_corr >= 0.4:
            strength = "moderate"
        elif abs_corr >= 0.2:
            strength = "weak"
        else:
            strength = "very weak"
            
        direction = "positive" if stats['correlation'] > 0 else "negative"
        significance = "significant" if stats['p_value'] < 0.05 else "not significant"
        
        print(f"  Interpretation: {strength} {direction} correlation ({significance})")
        print()

def main():
    """Main analysis function."""
    # Load data
    print("Loading telomere analysis data...")
    df = load_telomere_data('telomere_analysis.csv')
    print(f"Loaded {len(df)} samples with {len(df.columns)} features")
    
    # Calculate signature scores
    print("Calculating COSMIC SBS signature scores...")
    signatures = calculate_sbs_signatures(df)
    
    # Add signature scores to dataframe
    for sig_name, sig_scores in signatures.items():
        df[sig_name] = sig_scores
    
    # Calculate correlations with age
    print("Calculating correlations with age...")
    correlations = calculate_correlations_with_age(df, signatures)
    
    # Rank signatures
    sig_rankings = rank_signatures(correlations)
    
    # Print results
    print_detailed_results(correlations, sig_rankings)
    
    # Create visualizations
    print("Creating visualizations...")
    fig = create_visualization(df, signatures, correlations)
    
    # Save plots
    plt.savefig('cosmic_sbs_signatures_vs_age.png', dpi=300, bbox_inches='tight')
    print("Saved plot as 'cosmic_sbs_signatures_vs_age.png'")
    
    # Save results to CSV
    results_df = pd.DataFrame([
        {
            'Signature': sig_name.replace('_score', ''),
            'Spearman_r': stats['correlation'],
            'P_value': stats['p_value'],
            'N_samples': stats['n_samples'],
            'Abs_correlation': abs(stats['correlation'])
        }
        for sig_name, stats in correlations.items()
    ])
    results_df = results_df.sort_values('Abs_correlation', ascending=False)
    results_df.to_csv('cosmic_sbs_signatures_correlations.csv', index=False)
    print("Saved correlation results as 'cosmic_sbs_signatures_correlations.csv'")
    
    # Show plot
    plt.show()
    
    return df, signatures, correlations, sig_rankings

if __name__ == "__main__":
    df, signatures, correlations, sig_rankings = main()
