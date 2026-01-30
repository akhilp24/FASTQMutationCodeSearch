import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
from scipy.optimize import curve_fit
import sys
from pathlib import Path

# Add analysis directory to path for imports
analysis_dir = Path(__file__).parent.parent.parent / "analysis"
sys.path.append(str(analysis_dir))

# Import utility functions
sys.path.append(str(Path(__file__).parent.parent))
from utils import check_data_loaded

def show_statistical_analysis_page():
    """Display the statistical analysis page."""
    
    st.title("📊 Statistical Analysis")
    st.markdown("Comprehensive statistical analysis including correlations, composite scores, and curve fitting.")
    
    # Check if data is loaded
    if not check_data_loaded():
        return
    
    df = st.session_state.analysis_data
    
    # Analysis tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "🔗 Correlations", 
        "📈 Composite Scores", 
        "📊 Heatmaps", 
        "📉 Curve Fitting"
    ])
    
    with tab1:
        show_correlation_analysis(df)
    
    with tab2:
        show_composite_score_analysis(df)
    
    with tab3:
        show_heatmap_analysis(df)
    
    with tab4:
        show_curve_fitting_analysis(df)

def show_correlation_analysis(df):
    """Display Spearman correlation analysis."""
    
    st.subheader("🔗 Spearman Correlation Analysis")
    st.markdown("Calculate Spearman correlations between all numeric variables and age.")
    
    # Get numeric columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if 'Age' not in numeric_cols:
        st.error("Age column not found in numeric data.")
        return
    
    # Remove Age from the list to avoid self-correlation
    analysis_cols = [col for col in numeric_cols if col != 'Age']
    
    # Filter data
    df_clean = df.dropna(subset=['Age'])
    
    if len(df_clean) < 2:
        st.error("Not enough data points for correlation analysis.")
        return
    
    # Calculate correlations
    correlations = []
    for col in analysis_cols:
        col_data = df_clean[col].dropna()
        if len(col_data) >= 2:
            # Get matching age data
            age_data = df_clean.loc[col_data.index, 'Age']
            if len(age_data) >= 2:
                corr, p_val = stats.spearmanr(age_data, col_data)
                correlations.append({
                    'Variable': col,
                    'Spearman_r': corr,
                    'P_value': p_val,
                    'N_samples': len(age_data),
                    'Significant': 'Yes' if p_val < 0.05 else 'No'
                })
    
    if correlations:
        corr_df = pd.DataFrame(correlations)
        corr_df = corr_df.sort_values('Spearman_r', key=abs, ascending=False)
        
        # Display results
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.subheader("📊 Correlation Results")
            st.dataframe(corr_df, use_container_width=True)
            
            # Download button
            csv = corr_df.to_csv(index=False)
            st.download_button(
                label="📥 Download Correlation Results",
                data=csv,
                file_name="spearman_correlations.csv",
                mime="text/csv"
            )
        
        with col2:
            st.subheader("📈 Top Correlations")
            
            # Show top 10 correlations
            top_corr = corr_df.head(10)
            
            fig = px.bar(
                top_corr,
                x='Spearman_r',
                y='Variable',
                orientation='h',
                title="Top 10 Correlations with Age",
                color='Spearman_r',
                color_continuous_scale='RdBu_r'
            )
            fig.update_layout(height=400)
            st.plotly_chart(fig, use_container_width=True)
        
        # Interactive correlation plots
        st.subheader("🔍 Interactive Correlation Plots")
        
        selected_var = st.selectbox(
            "Select variable for detailed correlation plot:",
            options=corr_df['Variable'].tolist(),
            index=0
        )
        
        if selected_var:
            plot_data = df_clean[['Age', selected_var]].dropna()
            
            if len(plot_data) > 1:
                # Create scatter plot with trendline
                fig = px.scatter(
                    plot_data,
                    x='Age',
                    y=selected_var,
                    title=f"{selected_var} vs Age",
                    hover_data=['Age', selected_var]
                )
                
                # Add trendline
                z = np.polyfit(plot_data['Age'], plot_data[selected_var], 1)
                p = np.poly1d(z)
                x_trend = np.linspace(plot_data['Age'].min(), plot_data['Age'].max(), 100)
                y_trend = p(x_trend)
                
                fig.add_trace(go.Scatter(
                    x=x_trend,
                    y=y_trend,
                    mode='lines',
                    name='Trendline',
                    line=dict(color='red', dash='dash')
                ))
                
                # Get correlation info
                corr_info = corr_df[corr_df['Variable'] == selected_var].iloc[0]
                fig.add_annotation(
                    text=f"r = {corr_info['Spearman_r']:.3f}<br>p = {corr_info['P_value']:.3e}",
                    xref="paper", yref="paper",
                    x=0.05, y=0.95,
                    showarrow=False,
                    bgcolor="white",
                    bordercolor="black",
                    borderwidth=1
                )
                
                st.plotly_chart(fig, use_container_width=True)
    else:
        st.warning("No valid correlations could be calculated.")

def show_composite_score_analysis(df):
    """Display composite score analysis."""
    
    st.subheader("📈 Composite Score Analysis")
    st.markdown("Calculate and analyze composite scores for age prediction.")
    
    # Define top features for composite score (from the original code)
    top_features = [
        'total_mutations_per_1k_reads',
        'g_strand_mutations_A>C_a1_per_1k',
        'g_strand_mutations_T>G_t2_per_1k',
        'g_strand_mutations_T>G_t1_per_1k',
        'g_strand_mutations_G>A_g3_per_1k',
    ]
    
    # Check which features are available
    available_features = [f for f in top_features if f in df.columns]
    
    if not available_features:
        st.warning("No composite score features found in the data.")
        return
    
    # Calculate composite score
    df_clean = df.dropna(subset=['Age'] + available_features)
    
    if len(df_clean) < 2:
        st.error("Not enough data points for composite score analysis.")
        return
    
    # Calculate composite score
    df_clean['composite_score'] = df_clean[available_features].mean(axis=1)
    
    # Calculate correlation
    corr, p_val = stats.spearmanr(df_clean['composite_score'], df_clean['Age'])
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📊 Composite Score Results")
        
        # Display metrics
        st.metric("Spearman Correlation", f"{corr:.3f}")
        st.metric("P-value", f"{p_val:.3e}")
        st.metric("Sample Size", len(df_clean))
        
        # Feature importance
        st.subheader("🎯 Feature Contributions")
        feature_contributions = []
        for feature in available_features:
            feature_corr, _ = stats.spearmanr(df_clean[feature], df_clean['Age'])
            feature_contributions.append({
                'Feature': feature,
                'Correlation with Age': f"{feature_corr:.3f}",
                'Mean Value': f"{df_clean[feature].mean():.3f}"
            })
        
        contrib_df = pd.DataFrame(feature_contributions)
        st.dataframe(contrib_df, use_container_width=True)
    
    with col2:
        st.subheader("📈 Composite Score Plot")
        
        # Create scatter plot
        fig = px.scatter(
            df_clean,
            x='Age',
            y='composite_score',
            title="Composite Score vs Age",
            hover_data=['FileName', 'Age', 'composite_score']
        )
        
        # Add trendline
        z = np.polyfit(df_clean['Age'], df_clean['composite_score'], 1)
        p = np.poly1d(z)
        x_trend = np.linspace(df_clean['Age'].min(), df_clean['Age'].max(), 100)
        y_trend = p(x_trend)
        
        fig.add_trace(go.Scatter(
            x=x_trend,
            y=y_trend,
            mode='lines',
            name='Trendline',
            line=dict(color='red', dash='dash')
        ))
        
        # Add correlation annotation
        fig.add_annotation(
            text=f"r = {corr:.3f}<br>p = {p_val:.3e}",
            xref="paper", yref="paper",
            x=0.05, y=0.95,
            showarrow=False,
            bgcolor="white",
            bordercolor="black",
            borderwidth=1
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    # Show composite score distribution
    st.subheader("📊 Composite Score Distribution")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # Histogram
        fig = px.histogram(
            df_clean,
            x='composite_score',
            nbins=20,
            title="Composite Score Distribution"
        )
        st.plotly_chart(fig, use_container_width=True)
    
    with col2:
        # Box plot by age groups
        df_clean['Age_Group'] = pd.cut(
            df_clean['Age'], 
            bins=3, 
            labels=['Young', 'Middle', 'Old']
        )
        
        fig = px.box(
            df_clean,
            x='Age_Group',
            y='composite_score',
            title="Composite Score by Age Group"
        )
        st.plotly_chart(fig, use_container_width=True)

def show_heatmap_analysis(df):
    """Display correlation heatmap analysis."""
    
    st.subheader("📊 Correlation Heatmap Analysis")
    st.markdown("Visualize pairwise correlations between variables.")
    
    # Get mutation columns
    mutation_cols = [col for col in df.columns if 'per_1k' in col]
    
    # Add other important columns
    important_cols = ['Age', 'Telomere_Length', 'Total_Reads']
    for col in important_cols:
        if col in df.columns and col not in mutation_cols:
            mutation_cols.append(col)
    
    # Limit to avoid too large heatmaps
    if len(mutation_cols) > 20:
        st.warning(f"Too many variables ({len(mutation_cols)}). Showing first 20.")
        mutation_cols = mutation_cols[:20]
    
    # Calculate correlation matrix
    df_clean = df[mutation_cols].dropna()
    
    if len(df_clean) < 2:
        st.error("Not enough data points for heatmap analysis.")
        return
    
    corr_matrix = df_clean.corr()
    
    # Create heatmap
    fig = px.imshow(
        corr_matrix,
        text_auto=True,
        aspect="auto",
        title="Correlation Heatmap",
        color_continuous_scale='RdBu_r'
    )
    
    fig.update_layout(height=600)
    st.plotly_chart(fig, use_container_width=True)
    
    # Show strongest correlations
    st.subheader("🔍 Strongest Correlations")
    
    # Get upper triangle of correlation matrix
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    corr_pairs = []
    
    for i in range(len(corr_matrix.columns)):
        for j in range(i+1, len(corr_matrix.columns)):
            var1 = corr_matrix.columns[i]
            var2 = corr_matrix.columns[j]
            corr_val = corr_matrix.iloc[i, j]
            corr_pairs.append({
                'Variable 1': var1,
                'Variable 2': var2,
                'Correlation': corr_val,
                'Abs_Correlation': abs(corr_val)
            })
    
    corr_pairs_df = pd.DataFrame(corr_pairs)
    corr_pairs_df = corr_pairs_df.sort_values('Abs_Correlation', ascending=False)
    
    # Show top correlations
    st.dataframe(corr_pairs_df.head(15), use_container_width=True)

def show_curve_fitting_analysis(df):
    """Display curve fitting analysis."""
    
    st.subheader("📉 Curve Fitting Analysis")
    st.markdown("Fit different curve types to identify the best relationship between variables and age.")
    
    # Define curve fitting functions
    def linear_func(x, a, b):
        return a * x + b
    
    def exponential_func(x, a, b, c):
        return a * np.exp(b * x) + c
    
    def logarithmic_func(x, a, b):
        return a * np.log(x + 1) + b
    
    def polynomial_func(x, a, b, c, d):
        return a * x**3 + b * x**2 + c * x + d
    
    def power_func(x, a, b, c):
        return a * (x + 1)**b + c
    
    # Select variable for curve fitting
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    numeric_cols = [col for col in numeric_cols if col != 'Age']
    
    selected_var = st.selectbox(
        "Select variable for curve fitting:",
        options=numeric_cols,
        index=0
    )
    
    if selected_var:
        # Prepare data
        plot_data = df[['Age', selected_var]].dropna()
        
        if len(plot_data) < 4:
            st.error("Not enough data points for curve fitting (need at least 4).")
            return
        
        x_data = plot_data['Age'].values
        y_data = plot_data[selected_var].values
        
        # Try different curve types
        curve_types = [
            ('Linear', linear_func, [1, 1]),
            ('Exponential', exponential_func, [1, 0.01, 1]),
            ('Logarithmic', logarithmic_func, [1, 1]),
            ('Polynomial', polynomial_func, [0.01, 0.1, 1, 1]),
            ('Power', power_func, [1, 0.5, 1])
        ]
        
        results = []
        best_fit = None
        best_r_squared = -np.inf
        
        # Create subplot
        fig = make_subplots(rows=1, cols=1)
        
        # Plot original data
        fig.add_trace(go.Scatter(
            x=x_data,
            y=y_data,
            mode='markers',
            name='Data points',
            marker=dict(color='black', size=8)
        ))
        
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
                    'Curve_Type': name,
                    'R_squared': r_squared,
                    'Parameters': popt.tolist()
                })
                
                # Plot the fitted curve
                x_smooth = np.linspace(x_data.min(), x_data.max(), 100)
                y_smooth = func(x_smooth, *popt)
                
                fig.add_trace(go.Scatter(
                    x=x_smooth,
                    y=y_smooth,
                    mode='lines',
                    name=f'{name} (R² = {r_squared:.3f})',
                    line=dict(color=colors[i], dash='dash')
                ))
                
                # Track best fit
                if r_squared > best_r_squared:
                    best_r_squared = r_squared
                    best_fit = (name, func, popt)
                    
            except Exception as e:
                st.warning(f"Could not fit {name} curve: {e}")
        
        # Highlight best fit
        if best_fit:
            name, func, popt = best_fit
            x_smooth = np.linspace(x_data.min(), x_data.max(), 100)
            y_smooth = func(x_smooth, *popt)
            
            fig.add_trace(go.Scatter(
                x=x_smooth,
                y=y_smooth,
                mode='lines',
                name=f'Best fit: {name} (R² = {best_r_squared:.3f})',
                line=dict(color='red', width=3)
            ))
        
        fig.update_layout(
            title=f"Curve Fitting: {selected_var} vs Age",
            xaxis_title="Age (years)",
            yaxis_title=selected_var,
            height=500
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Display results table
        if results:
            results_df = pd.DataFrame(results)
            results_df = results_df.sort_values('R_squared', ascending=False)
            
            st.subheader("📊 Curve Fitting Results")
            st.dataframe(results_df, use_container_width=True)
            
            # Show best fit summary
            best_result = results_df.iloc[0]
            st.success(f"🏆 Best fit: {best_result['Curve_Type']} with R² = {best_result['R_squared']:.4f}")
            
            # Download results
            csv = results_df.to_csv(index=False)
            st.download_button(
                label="📥 Download Curve Fitting Results",
                data=csv,
                file_name=f"curve_fitting_{selected_var}.csv",
                mime="text/csv"
            )

# Call the function when the page is loaded
show_statistical_analysis_page()
