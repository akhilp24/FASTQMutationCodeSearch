import streamlit as st
import pandas as pd
import numpy as np
import os
import zipfile
import io
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
import sys

# Add analysis directory to path for imports
analysis_dir = Path(__file__).parent.parent.parent / "analysis"
sys.path.append(str(analysis_dir))

# Import utility functions
sys.path.append(str(Path(__file__).parent.parent))
from utils import check_data_loaded

def show_export_results_page():
    """Display the export results page."""
    
    st.title("💾 Export Results")
    st.markdown("Download your analysis results, plots, and data in various formats.")
    
    # Check if data is loaded
    if not check_data_loaded():
        return
    
    df = st.session_state.analysis_data
    
    # Export tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "📊 Data Export", 
        "📈 Plot Export", 
        "📋 Report Generation", 
        "🗂️ Batch Download"
    ])
    
    with tab1:
        show_data_export(df)
    
    with tab2:
        show_plot_export(df)
    
    with tab3:
        show_report_generation(df)
    
    with tab4:
        show_batch_download(df)

def show_data_export(df):
    """Display data export options."""
    
    st.subheader("📊 Data Export")
    st.markdown("Download your analysis data in various formats.")
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📋 Raw Data")
        
        # Full dataset
        csv_data = df.to_csv(index=False)
        st.download_button(
            label="📥 Download Full Dataset (CSV)",
            data=csv_data,
            file_name="telomere_analysis_full.csv",
            mime="text/csv"
        )
        
        # Filtered data options
        st.subheader("🔍 Filtered Data")
        
        # Age range filter
        age_min, age_max = st.slider(
            "Age Range",
            min_value=int(df['Age'].min()),
            max_value=int(df['Age'].max()),
            value=(int(df['Age'].min()), int(df['Age'].max()))
        )
        
        # Sample selection
        selected_samples = st.multiselect(
            "Select Samples",
            options=df['FileName'].tolist(),
            default=df['FileName'].tolist()
        )
        
        # Filter data
        filtered_df = df[
            (df['Age'] >= age_min) & 
            (df['Age'] <= age_max) & 
            (df['FileName'].isin(selected_samples))
        ]
        
        if len(filtered_df) != len(df):
            st.info(f"Filtered to {len(filtered_df)} samples (from {len(df)} total)")
            
            filtered_csv = filtered_df.to_csv(index=False)
            st.download_button(
                label="📥 Download Filtered Dataset (CSV)",
                data=filtered_csv,
                file_name="telomere_analysis_filtered.csv",
                mime="text/csv"
            )
    
    with col2:
        st.subheader("📊 Summary Statistics")
        
        # Generate summary statistics
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        summary_stats = df[numeric_cols].describe()
        
        # Add additional statistics
        summary_stats.loc['skewness'] = df[numeric_cols].skew()
        summary_stats.loc['kurtosis'] = df[numeric_cols].kurtosis()
        
        st.dataframe(summary_stats, use_container_width=True)
        
        # Download summary
        summary_csv = summary_stats.to_csv()
        st.download_button(
            label="📥 Download Summary Statistics (CSV)",
            data=summary_csv,
            file_name="summary_statistics.csv",
            mime="text/csv"
        )
        
        # Correlation matrix
        st.subheader("🔗 Correlation Matrix")
        
        if len(numeric_cols) > 1:
            corr_matrix = df[numeric_cols].corr()
            
            # Show heatmap
            fig = px.imshow(
                corr_matrix,
                text_auto=True,
                aspect="auto",
                title="Correlation Matrix"
            )
            fig.update_layout(height=400)
            st.plotly_chart(fig, use_container_width=True)
            
            # Download correlation matrix
            corr_csv = corr_matrix.to_csv()
            st.download_button(
                label="📥 Download Correlation Matrix (CSV)",
                data=corr_csv,
                file_name="correlation_matrix.csv",
                mime="text/csv"
            )

def show_plot_export(df):
    """Display plot export options."""
    
    st.subheader("📈 Plot Export")
    st.markdown("Generate and download publication-ready plots.")
    
    # Plot type selection
    plot_type = st.selectbox(
        "Select plot type:",
        [
            "Scatter Plot (Age vs Variable)",
            "Correlation Heatmap",
            "Distribution Histograms",
            "Box Plots by Age Groups",
            "Trendline Analysis",
            "Mutational Signatures"
        ]
    )
    
    if plot_type == "Scatter Plot (Age vs Variable)":
        export_scatter_plot(df)
    elif plot_type == "Correlation Heatmap":
        export_correlation_heatmap(df)
    elif plot_type == "Distribution Histograms":
        export_distribution_plots(df)
    elif plot_type == "Box Plots by Age Groups":
        export_box_plots(df)
    elif plot_type == "Trendline Analysis":
        export_trendline_plots(df)
    elif plot_type == "Mutational Signatures":
        export_signature_plots(df)

def export_scatter_plot(df):
    """Export scatter plot."""
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        # Variable selection
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        y_var = st.selectbox("Y-axis variable:", numeric_cols)
        
        # Plot options
        show_trendline = st.checkbox("Show trendline", value=True)
        show_correlation = st.checkbox("Show correlation info", value=True)
        
        # Color options
        color_by = st.selectbox("Color by:", ['None'] + numeric_cols)
        
        # Export format
        export_format = st.selectbox("Export format:", ["PNG", "PDF", "SVG"])
    
    with col2:
        # Create plot
        fig = px.scatter(
            df,
            x='Age',
            y=y_var,
            color=color_by if color_by != 'None' else None,
            title=f"{y_var} vs Age",
            labels={'Age': 'Age (years)', y_var: y_var}
        )
        
        if show_trendline:
            # Add trendline
            plot_data = df[['Age', y_var]].dropna()
            if len(plot_data) > 1:
                z = np.polyfit(plot_data['Age'], plot_data[y_var], 1)
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
        
        if show_correlation:
            plot_data = df[['Age', y_var]].dropna()
            if len(plot_data) > 1:
                corr, p_val = stats.spearmanr(plot_data['Age'], plot_data[y_var])
                fig.add_annotation(
                    text=f"r = {corr:.3f}<br>p = {p_val:.3e}",
                    xref="paper", yref="paper",
                    x=0.05, y=0.95,
                    showarrow=False,
                    bgcolor="white",
                    bordercolor="black",
                    borderwidth=1
                )
        
        fig.update_layout(height=500)
        st.plotly_chart(fig, use_container_width=True)
        
        # Export button
        if export_format == "PNG":
            img_bytes = fig.to_image(format="png", width=800, height=600, scale=2)
            st.download_button(
                label="📥 Download Plot (PNG)",
                data=img_bytes,
                file_name=f"scatter_plot_{y_var}.png",
                mime="image/png"
            )
        elif export_format == "PDF":
            img_bytes = fig.to_image(format="pdf", width=800, height=600)
            st.download_button(
                label="📥 Download Plot (PDF)",
                data=img_bytes,
                file_name=f"scatter_plot_{y_var}.pdf",
                mime="application/pdf"
            )
        elif export_format == "SVG":
            img_bytes = fig.to_image(format="svg", width=800, height=600)
            st.download_button(
                label="📥 Download Plot (SVG)",
                data=img_bytes,
                file_name=f"scatter_plot_{y_var}.svg",
                mime="image/svg+xml"
            )

def export_correlation_heatmap(df):
    """Export correlation heatmap."""
    
    # Get numeric columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    
    # Limit columns for readability
    if len(numeric_cols) > 20:
        st.warning(f"Too many variables ({len(numeric_cols)}). Showing first 20.")
        numeric_cols = numeric_cols[:20]
    
    # Calculate correlation matrix
    corr_matrix = df[numeric_cols].corr()
    
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
    
    # Export options
    col1, col2, col3 = st.columns(3)
    
    with col1:
        img_bytes = fig.to_image(format="png", width=1000, height=800, scale=2)
        st.download_button(
            label="📥 Download PNG",
            data=img_bytes,
            file_name="correlation_heatmap.png",
            mime="image/png"
        )
    
    with col2:
        img_bytes = fig.to_image(format="pdf", width=1000, height=800)
        st.download_button(
            label="📥 Download PDF",
            data=img_bytes,
            file_name="correlation_heatmap.pdf",
            mime="application/pdf"
        )
    
    with col3:
        # Also export correlation matrix as CSV
        corr_csv = corr_matrix.to_csv()
        st.download_button(
            label="📥 Download CSV",
            data=corr_csv,
            file_name="correlation_matrix.csv",
            mime="text/csv"
        )

def export_distribution_plots(df):
    """Export distribution plots."""
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        # Variable selection
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        selected_vars = st.multiselect(
            "Select variables for distribution plots:",
            options=numeric_cols,
            default=numeric_cols[:4] if len(numeric_cols) >= 4 else numeric_cols
        )
        
        plot_type = st.selectbox("Plot type:", ["Histogram", "Box Plot", "Violin Plot"])
    
    with col2:
        if selected_vars:
            # Create subplots
            n_vars = len(selected_vars)
            n_cols = min(2, n_vars)
            n_rows = (n_vars + 1) // 2
            
            fig = make_subplots(
                rows=n_rows, cols=n_cols,
                subplot_titles=selected_vars,
                vertical_spacing=0.1
            )
            
            for i, var in enumerate(selected_vars):
                row = (i // n_cols) + 1
                col = (i % n_cols) + 1
                
                if plot_type == "Histogram":
                    fig.add_trace(
                        go.Histogram(x=df[var], name=var, showlegend=False),
                        row=row, col=col
                    )
                elif plot_type == "Box Plot":
                    fig.add_trace(
                        go.Box(y=df[var], name=var, showlegend=False),
                        row=row, col=col
                    )
                elif plot_type == "Violin Plot":
                    fig.add_trace(
                        go.Violin(y=df[var], name=var, showlegend=False),
                        row=row, col=col
                    )
            
            fig.update_layout(
                height=400 * n_rows,
                title_text=f"Distribution Plots - {plot_type}",
                showlegend=False
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Export button
            img_bytes = fig.to_image(format="png", width=1000, height=400*n_rows, scale=2)
            st.download_button(
                label="📥 Download Distribution Plots (PNG)",
                data=img_bytes,
                file_name=f"distribution_plots_{plot_type.lower()}.png",
                mime="image/png"
            )

def export_box_plots(df):
    """Export box plots by age groups."""
    
    # Create age groups
    df_copy = df.copy()
    df_copy['Age_Group'] = pd.cut(
        df_copy['Age'], 
        bins=3, 
        labels=['Young', 'Middle', 'Old']
    )
    
    # Variable selection
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    numeric_cols = [col for col in numeric_cols if col != 'Age']
    
    selected_vars = st.multiselect(
        "Select variables for box plots:",
        options=numeric_cols,
        default=numeric_cols[:4] if len(numeric_cols) >= 4 else numeric_cols
    )
    
    if selected_vars:
        # Create subplots
        n_vars = len(selected_vars)
        n_cols = min(2, n_vars)
        n_rows = (n_vars + 1) // 2
        
        fig = make_subplots(
            rows=n_rows, cols=n_cols,
            subplot_titles=selected_vars,
            vertical_spacing=0.1
        )
        
        for i, var in enumerate(selected_vars):
            row = (i // n_cols) + 1
            col = (i % n_cols) + 1
            
            for age_group in ['Young', 'Middle', 'Old']:
                group_data = df_copy[df_copy['Age_Group'] == age_group][var].dropna()
                if len(group_data) > 0:
                    fig.add_trace(
                        go.Box(
                            y=group_data,
                            name=age_group,
                            legendgroup=age_group,
                            showlegend=(i == 0)  # Only show legend for first subplot
                        ),
                        row=row, col=col
                    )
        
        fig.update_layout(
            height=400 * n_rows,
            title_text="Box Plots by Age Groups",
            showlegend=True
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Export button
        img_bytes = fig.to_image(format="png", width=1000, height=400*n_rows, scale=2)
        st.download_button(
            label="📥 Download Box Plots (PNG)",
            data=img_bytes,
            file_name="box_plots_age_groups.png",
            mime="image/png"
        )

def export_trendline_plots(df):
    """Export trendline analysis plots."""
    
    # Get mutation rate columns
    mutation_cols = [col for col in df.columns if 'per_1k' in col and 'total' in col.lower()]
    
    if not mutation_cols:
        st.warning("No mutation rate columns found for trendline analysis.")
        return
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[col.replace('_', ' ').title() for col in mutation_cols[:4]],
        vertical_spacing=0.1
    )
    
    colors = ['blue', 'green', 'red', 'orange']
    
    for i, col in enumerate(mutation_cols[:4]):
        row = (i // 2) + 1
        col_idx = (i % 2) + 1
        
        # Filter out NaN values
        plot_data = df[['Age', col]].dropna()
        
        if len(plot_data) > 1:
            # Add scatter plot
            fig.add_trace(
                go.Scatter(
                    x=plot_data['Age'],
                    y=plot_data[col],
                    mode='markers',
                    name=col,
                    marker=dict(color=colors[i], size=8),
                    showlegend=False
                ),
                row=row, col=col_idx
            )
            
            # Add trendline
            z = np.polyfit(plot_data['Age'], plot_data[col], 1)
            p = np.poly1d(z)
            x_trend = np.linspace(plot_data['Age'].min(), plot_data['Age'].max(), 100)
            y_trend = p(x_trend)
            
            fig.add_trace(
                go.Scatter(
                    x=x_trend,
                    y=y_trend,
                    mode='lines',
                    name=f'{col} trend',
                    line=dict(color=colors[i], dash='dash'),
                    showlegend=False
                ),
                row=row, col=col_idx
            )
            
            # Calculate R-squared
            y_pred = p(plot_data['Age'])
            ss_res = np.sum((plot_data[col] - y_pred) ** 2)
            ss_tot = np.sum((plot_data[col] - np.mean(plot_data[col])) ** 2)
            r_squared = 1 - (ss_res / ss_tot)
            
            # Add R-squared annotation
            fig.add_annotation(
                text=f"R² = {r_squared:.3f}",
                xref=f"x{i+1}", yref=f"y{i+1}",
                x=0.05, y=0.95,
                showarrow=False,
                bgcolor="white",
                bordercolor="black",
                borderwidth=1
            )
            
            # Update axis labels
            fig.update_xaxes(title_text="Age (years)", row=row, col=col_idx)
            fig.update_yaxes(title_text="Mutations per 1k", row=row, col=col_idx)
    
    fig.update_layout(
        height=600,
        title_text="Trendline Analysis: Mutation Rates vs Age",
        showlegend=False
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Export button
    img_bytes = fig.to_image(format="png", width=1200, height=800, scale=2)
    st.download_button(
        label="📥 Download Trendline Plots (PNG)",
        data=img_bytes,
        file_name="trendline_analysis.png",
        mime="image/png"
    )

def export_signature_plots(df):
    """Export mutational signature plots."""
    
    # Calculate signatures (simplified version)
    signatures = {}
    
    # SBS1: C>T mutations
    ct_mutations = 0
    if 'c_strand_C>T_per_1k' in df.columns:
        ct_mutations += df['c_strand_C>T_per_1k']
    if 'g_strand_G>A_per_1k' in df.columns:
        ct_mutations += df['g_strand_G>A_per_1k']
    signatures['SBS1_CT_score'] = ct_mutations
    
    # SBS4: C>A mutations
    ca_mutations = 0
    if 'c_strand_C>A_per_1k' in df.columns:
        ca_mutations += df['c_strand_C>A_per_1k']
    if 'g_strand_G>T_per_1k' in df.columns:
        ca_mutations += df['g_strand_G>T_per_1k']
    signatures['SBS4_CA_score'] = ca_mutations
    
    # Only use samples with age data
    age_mask = df['Age'].notna()
    age_data = df[age_mask]['Age']
    
    if len(age_data) < 2:
        st.error("Not enough age data for signature plots.")
        return
    
    # Create subplots
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=['SBS1: C>T (Clock-like)', 'SBS4: C>A (Tobacco smoking)'],
        vertical_spacing=0.1
    )
    
    sig_names = ['SBS1_CT_score', 'SBS4_CA_score']
    colors = ['blue', 'green']
    
    for i, sig_name in enumerate(sig_names):
        if sig_name in signatures:
            sig_data = signatures[sig_name][age_mask]
            
            # Create scatter plot
            fig.add_trace(
                go.Scatter(
                    x=age_data,
                    y=sig_data,
                    mode='markers',
                    name=sig_name,
                    marker=dict(color=colors[i], size=8),
                    showlegend=False
                ),
                row=1, col=i+1
            )
            
            # Add trend line
            if len(sig_data) > 1:
                z = np.polyfit(age_data, sig_data, 1)
                p = np.poly1d(z)
                x_trend = np.linspace(age_data.min(), age_data.max(), 100)
                y_trend = p(x_trend)
                
                fig.add_trace(
                    go.Scatter(
                        x=x_trend,
                        y=y_trend,
                        mode='lines',
                        name=f'{sig_name} trend',
                        line=dict(color=colors[i], dash='dash'),
                        showlegend=False
                    ),
                    row=1, col=i+1
                )
                
                # Calculate correlation
                corr, p_val = stats.spearmanr(age_data, sig_data)
                
                # Add correlation info
                fig.add_annotation(
                    text=f"r = {corr:.3f}<br>p = {p_val:.3e}",
                    xref=f"x{i+1}", yref=f"y{i+1}",
                    x=0.05, y=0.95,
                    showarrow=False,
                    bgcolor="white",
                    bordercolor="black",
                    borderwidth=1
                )
            
            # Update axis labels
            fig.update_xaxes(title_text="Age (years)", row=1, col=i+1)
            fig.update_yaxes(title_text="Mutations per 1k bases", row=1, col=i+1)
    
    fig.update_layout(
        height=400,
        title_text="COSMIC SBS Mutational Signatures vs Age",
        showlegend=False
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Export button
    img_bytes = fig.to_image(format="png", width=1200, height=600, scale=2)
    st.download_button(
        label="📥 Download Signature Plots (PNG)",
        data=img_bytes,
        file_name="mutational_signatures.png",
        mime="image/png"
    )

def show_report_generation(df):
    """Display report generation options."""
    
    st.subheader("📋 Report Generation")
    st.markdown("Generate comprehensive analysis reports.")
    
    # Report options
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📊 Report Contents")
        
        include_summary = st.checkbox("Summary Statistics", value=True)
        include_correlations = st.checkbox("Correlation Analysis", value=True)
        include_plots = st.checkbox("Key Plots", value=True)
        include_signatures = st.checkbox("Mutational Signatures", value=True)
        include_metadata = st.checkbox("Sample Metadata", value=True)
    
    with col2:
        st.subheader("📈 Report Format")
        
        report_format = st.selectbox("Report format:", ["HTML", "PDF", "Markdown"])
        include_data = st.checkbox("Include raw data", value=False)
        plot_resolution = st.selectbox("Plot resolution:", ["High (300 DPI)", "Medium (150 DPI)", "Low (72 DPI)"])
    
    # Generate report button
    if st.button("📋 Generate Report", type="primary"):
        generate_analysis_report(df, {
            'include_summary': include_summary,
            'include_correlations': include_correlations,
            'include_plots': include_plots,
            'include_signatures': include_signatures,
            'include_metadata': include_metadata,
            'include_data': include_data,
            'format': report_format,
            'resolution': plot_resolution
        })

def generate_analysis_report(df, options):
    """Generate analysis report."""
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    try:
        status_text.text("📋 Generating report...")
        progress_bar.progress(0.2)
        
        # Create report content
        report_content = []
        
        # Header
        report_content.append("# Telomere Analysis Report")
        report_content.append(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_content.append(f"Total samples: {len(df)}")
        report_content.append("")
        
        progress_bar.progress(0.4)
        
        # Summary statistics
        if options['include_summary']:
            report_content.append("## Summary Statistics")
            numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
            summary_stats = df[numeric_cols].describe()
            report_content.append(summary_stats.to_string())
            report_content.append("")
        
        progress_bar.progress(0.6)
        
        # Correlation analysis
        if options['include_correlations']:
            report_content.append("## Correlation Analysis")
            if 'Age' in df.columns:
                age_correlations = []
                for col in df.select_dtypes(include=[np.number]).columns:
                    if col != 'Age':
                        corr_data = df[['Age', col]].dropna()
                        if len(corr_data) > 1:
                            corr, p_val = stats.spearmanr(corr_data['Age'], corr_data[col])
                            age_correlations.append(f"{col}: r = {corr:.3f}, p = {p_val:.3e}")
                
                report_content.extend(age_correlations[:10])  # Top 10 correlations
                report_content.append("")
        
        progress_bar.progress(0.8)
        
        # Sample metadata
        if options['include_metadata']:
            report_content.append("## Sample Information")
            report_content.append(f"Age range: {df['Age'].min():.1f} - {df['Age'].max():.1f} years")
            report_content.append(f"Mean age: {df['Age'].mean():.1f} years")
            if 'Telomere_Length' in df.columns:
                report_content.append(f"Telomere length range: {df['Telomere_Length'].min():.0f} - {df['Telomere_Length'].max():.0f} bp")
                report_content.append(f"Mean telomere length: {df['Telomere_Length'].mean():.0f} bp")
            report_content.append("")
        
        progress_bar.progress(1.0)
        status_text.text("✅ Report generated!")
        
        # Convert to requested format
        if options['format'] == "Markdown":
            report_text = "\n".join(report_content)
            st.download_button(
                label="📥 Download Report (Markdown)",
                data=report_text,
                file_name="telomere_analysis_report.md",
                mime="text/markdown"
            )
        elif options['format'] == "HTML":
            # Convert markdown to HTML (simplified)
            html_content = "<html><body>"
            for line in report_content:
                if line.startswith("# "):
                    html_content += f"<h1>{line[2:]}</h1>"
                elif line.startswith("## "):
                    html_content += f"<h2>{line[3:]}</h2>"
                elif line.startswith("### "):
                    html_content += f"<h3>{line[4:]}</h3>"
                elif line == "":
                    html_content += "<br>"
                else:
                    html_content += f"<p>{line}</p>"
            html_content += "</body></html>"
            
            st.download_button(
                label="📥 Download Report (HTML)",
                data=html_content,
                file_name="telomere_analysis_report.html",
                mime="text/html"
            )
        
        # Display preview
        st.subheader("📋 Report Preview")
        st.text("\n".join(report_content[:50]))  # Show first 50 lines
        
    except Exception as e:
        st.error(f"Error generating report: {e}")
        progress_bar.progress(0)
        status_text.text("❌ Report generation failed")

def show_batch_download(df):
    """Display batch download options."""
    
    st.subheader("🗂️ Batch Download")
    st.markdown("Download all analysis results in a single zip file.")
    
    # Download options
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📊 Data Files")
        
        include_full_data = st.checkbox("Full dataset (CSV)", value=True)
        include_filtered_data = st.checkbox("Filtered dataset (CSV)", value=False)
        include_summary_stats = st.checkbox("Summary statistics (CSV)", value=True)
        include_correlation_matrix = st.checkbox("Correlation matrix (CSV)", value=True)
    
    with col2:
        st.subheader("📈 Plot Files")
        
        include_scatter_plots = st.checkbox("Scatter plots (PNG)", value=True)
        include_heatmaps = st.checkbox("Correlation heatmaps (PNG)", value=True)
        include_distributions = st.checkbox("Distribution plots (PNG)", value=True)
        include_signature_plots = st.checkbox("Signature plots (PNG)", value=True)
    
    # Generate batch download
    if st.button("🗂️ Create Batch Download", type="primary"):
        create_batch_download(df, {
            'include_full_data': include_full_data,
            'include_filtered_data': include_filtered_data,
            'include_summary_stats': include_summary_stats,
            'include_correlation_matrix': include_correlation_matrix,
            'include_scatter_plots': include_scatter_plots,
            'include_heatmaps': include_heatmaps,
            'include_distributions': include_distributions,
            'include_signature_plots': include_signature_plots
        })

def create_batch_download(df, options):
    """Create batch download zip file."""
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    try:
        status_text.text("🗂️ Creating batch download...")
        progress_bar.progress(0.1)
        
        # Create zip file in memory
        zip_buffer = io.BytesIO()
        
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            file_count = 0
            
            # Add data files
            if options['include_full_data']:
                csv_data = df.to_csv(index=False)
                zip_file.writestr("data/full_dataset.csv", csv_data)
                file_count += 1
            
            if options['include_summary_stats']:
                numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
                summary_stats = df[numeric_cols].describe()
                summary_csv = summary_stats.to_csv()
                zip_file.writestr("data/summary_statistics.csv", summary_csv)
                file_count += 1
            
            if options['include_correlation_matrix']:
                numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
                corr_matrix = df[numeric_cols].corr()
                corr_csv = corr_matrix.to_csv()
                zip_file.writestr("data/correlation_matrix.csv", corr_csv)
                file_count += 1
            
            progress_bar.progress(0.5)
            
            # Add plot files (simplified - would need actual plot generation)
            if options['include_scatter_plots']:
                # Create a simple scatter plot
                fig = px.scatter(df, x='Age', y='Telomere_Length', title='Age vs Telomere Length')
                img_bytes = fig.to_image(format="png", width=800, height=600, scale=2)
                zip_file.writestr("plots/scatter_plot_age_vs_telomere.png", img_bytes)
                file_count += 1
            
            progress_bar.progress(0.8)
            
            # Add README
            readme_content = f"""
# Telomere Analysis Results

This zip file contains the results of your telomere analysis.

## Contents:
- data/: CSV files with analysis results
- plots/: PNG files with visualizations

## Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
## Total samples: {len(df)}
## Files included: {file_count}

## File Descriptions:
- full_dataset.csv: Complete analysis dataset
- summary_statistics.csv: Statistical summary of all variables
- correlation_matrix.csv: Correlation matrix between all numeric variables
- scatter_plot_age_vs_telomere.png: Scatter plot of age vs telomere length

For questions about this analysis, please refer to the Telomere Analysis Dashboard.
"""
            zip_file.writestr("README.txt", readme_content)
        
        progress_bar.progress(1.0)
        status_text.text("✅ Batch download ready!")
        
        # Download button
        zip_buffer.seek(0)
        st.download_button(
            label="📥 Download All Results (ZIP)",
            data=zip_buffer.getvalue(),
            file_name=f"telomere_analysis_results_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.zip",
            mime="application/zip"
        )
        
    except Exception as e:
        st.error(f"Error creating batch download: {e}")
        progress_bar.progress(0)
        status_text.text("❌ Batch download failed")

# Call the function when the page is loaded
show_export_results_page()
