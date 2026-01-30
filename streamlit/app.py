import streamlit as st
import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

# Add the analysis directory to Python path to import analysis modules
analysis_dir = Path(__file__).parent.parent / "analysis"
sys.path.append(str(analysis_dir))

# Page configuration
st.set_page_config(
    page_title="Telomere Mutational Quantification Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
    }
    .sidebar .sidebar-content {
        background-color: #f8f9fa;
    }
    .stButton > button {
        width: 100%;
        border-radius: 0.5rem;
    }
</style>
""", unsafe_allow_html=True)

# Initialize session state
if 'data_loaded' not in st.session_state:
    st.session_state.data_loaded = False
if 'analysis_data' not in st.session_state:
    st.session_state.analysis_data = None
if 'csv_generated' not in st.session_state:
    st.session_state.csv_generated = False

# Main header
st.markdown('<h2 class="main-header">🧬 Telomere Mutational Quantification Analysis Dashboard</h2>', unsafe_allow_html=True)


col1, col2 = st.columns([2, 1])

with col1:
    st.markdown("""
    ### About This Project
    
    This interactive dashboard provides comprehensive analysis tools for telomere sequencing data, 
    including:
    
    - **Data Processing**: Upload and process FASTQ/FASTA files to generate analysis-ready CSV files
    - **Visualizations**: Interactive plots for mutation patterns, correlations, and trends
    - **Statistical Analysis**: Spearman correlations, composite scores, and curve fitting
    - **Mutational Signatures**: COSMIC SBS signature analysis and correlation with age
    - **Export Tools**: Download plots, results, and analysis reports
    
    ### Getting Started
    
    1. **Upload Data**: Go to the Data Processing page to upload your FASTQ/FASTA files
    2. **Generate Analysis**: Process your files to create the analysis CSV
    3. **Explore Results**: Use the visualization and analysis pages to explore your data
    4. **Export Findings**: Download plots and results for your research
    
    ### Analysis Methods
    
    The dashboard implements sophisticated telomere analysis methods including:
    - Strand-specific mutation counting and normalization
    - Composite score calculations for age correlation
    - COSMIC mutational signature analysis
    - Statistical correlation analysis with multiple curve fitting models
    """)

with col2:
    st.markdown("### Quick Stats")
    
    # Check if data is available
    analysis_dir = Path(__file__).parent.parent / "analysis"
    csv_files = list(analysis_dir.glob('telomere_analysis*.csv'))
    if csv_files:
        try:
            df = pd.read_csv(csv_files[0])
            st.markdown(f"""
            <div class="metric-card">
                <h4>📊 Current Dataset</h4>
                <p><strong>Samples:</strong> {len(df)}</p>
                <p><strong>Features:</strong> {len(df.columns)}</p>
                <p><strong>Age Range:</strong> {df['Age'].min():.0f} - {df['Age'].max():.0f} years</p>
            </div>
            """, unsafe_allow_html=True)
            
            st.session_state.data_loaded = True
            st.session_state.analysis_data = df
            
        except Exception as e:
            st.error(f"Error loading data: {e}")
    else:
        st.markdown("""
        <div class="metric-card">
            <h4>📊 No Data Loaded</h4>
            <p>Upload and process your data to see statistics here.</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Show available analysis files
    st.markdown("### Available Analysis Files")
    analysis_files = [
        "generate_csv.py - CSV generation from sequencing files",
        "plotting.py - Comprehensive plotting functions", 
        "scatterplot.py - Scatter plot visualizations",
        "trendline.py - Trendline and correlation analysis",
        "histogram.py - Histogram and distribution analysis",
        "signatures.py - COSMIC SBS mutational signatures"
    ]
    
    for file_desc in analysis_files:
        st.markdown(f"• {file_desc}")

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666;'>
    <p>Telomere Analysis Dashboard | Built by Akhil Peddikuppa</p>
</div>
""", unsafe_allow_html=True)