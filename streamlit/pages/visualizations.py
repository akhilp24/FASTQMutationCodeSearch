import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
from pathlib import Path

# Add analysis directory to path for imports
analysis_dir = Path(__file__).parent.parent.parent / "analysis"
sys.path.append(str(analysis_dir))

# Import utility functions
sys.path.append(str(Path(__file__).parent.parent))
from utils import check_data_loaded

def show_visualizations_page():
    """Display the interactive visualizations page."""
    
    st.title("📈 Interactive Visualizations")
    st.markdown("Explore your telomere analysis data through interactive plots and charts.")
    
    # Check if data is loaded
    if not check_data_loaded():
        return
    
    df = st.session_state.analysis_data
    
    # Sidebar controls
    st.sidebar.header("🎛️ Plot Controls")
    
    # Age range filter
    age_min, age_max = st.sidebar.slider(
        "Age Range",
        min_value=int(df['Age'].min()),
        max_value=int(df['Age'].max()),
        value=(int(df['Age'].min()), int(df['Age'].max())),
        help="Filter data by age range"
    )
    
    # Sample selection
    selected_samples = st.sidebar.multiselect(
        "Select Samples",
        options=df['FileName'].tolist(),
        default=df['FileName'].tolist(),
        help="Choose specific samples to include in plots"
    )
    
    # Filter data based on selections
    filtered_df = df[
        (df['Age'] >= age_min) & 
        (df['Age'] <= age_max) & 
        (df['FileName'].isin(selected_samples))
    ]
    
    st.info(f"📊 Showing {len(filtered_df)} samples (filtered from {len(df)} total)")
    
    # Visualization tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "🔍 Scatter Plots", 
        "📊 Trendlines & Correlations", 
        "📈 Histograms", 
        "🎨 Custom Plots"
    ])
    
    with tab1:
        show_scatter_plots(filtered_df)
    
    with tab2:
        show_trendline_plots(filtered_df)
    
    with tab3:
        show_histogram_plots(filtered_df)
    
    with tab4:
        show_custom_plots(filtered_df)

def show_scatter_plots(df):
    """Display interactive scatter plots."""
    
    st.subheader("🔍 Scatter Plot Analysis")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        # Plot configuration
        x_var = st.selectbox(
            "X-axis variable:",
            options=['Age', 'Telomere_Length', 'Total_Reads'],
            index=0
        )
        
        # Get mutation columns for Y-axis
        mutation_cols = [col for col in df.columns if 'per_1k' in col]
        y_var = st.selectbox(
            "Y-axis variable:",
            options=mutation_cols,
            index=0
        )
        
        color_var = st.selectbox(
            "Color by:",
            options=['Age', 'Telomere_Length', 'Total_Reads', 'None'],
            index=0
        )
        
        size_var = st.selectbox(
            "Size by:",
            options=['Total_Reads', 'Age', 'Telomere_Length', 'None'],
            index=0
        )
        
        show_trendline = st.checkbox("Show trendline", value=True)
        show_labels = st.checkbox("Show sample labels", value=False)
    
    with col2:
        # Create interactive plot
        fig = px.scatter(
            df,
            x=x_var,
            y=y_var,
            color=color_var if color_var != 'None' else None,
            size=size_var if size_var != 'None' else None,
            hover_data=['FileName', 'Age', 'Telomere_Length'],
            title=f"{y_var} vs {x_var}",
            labels={x_var: x_var, y_var: y_var}
        )
        
        if show_trendline:
            # Add trendline
            z = np.polyfit(df[x_var].dropna(), df[y_var].dropna(), 1)
            p = np.poly1d(z)
            x_trend = np.linspace(df[x_var].min(), df[x_var].max(), 100)
            y_trend = p(x_trend)
            
            fig.add_trace(go.Scatter(
                x=x_trend,
                y=y_trend,
                mode='lines',
                name='Trendline',
                line=dict(color='red', dash='dash')
            ))
        
        if show_labels:
            fig.update_traces(
                text=df['FileName'],
                textposition="top center"
            )
        
        fig.update_layout(
            height=500,
            showlegend=True
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Calculate correlation
        if len(df) > 1:
            corr = df[x_var].corr(df[y_var])
            st.metric("Correlation Coefficient", f"{corr:.3f}")

def show_trendline_plots(df):
    """Display trendline analysis plots."""
    
    st.subheader("📊 Trendline & Correlation Analysis")
    
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
            
            # Update subplot title with R²
            fig.update_xaxes(title_text="Age (years)", row=row, col=col_idx)
            fig.update_yaxes(title_text="Mutations per 1k", row=row, col=col_idx)
    
    fig.update_layout(
        height=600,
        title_text="Mutation Rates vs Age with Trendlines",
        showlegend=False
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Correlation summary table
    st.subheader("📊 Correlation Summary")
    correlations = []
    for col in mutation_cols:
        plot_data = df[['Age', col]].dropna()
        if len(plot_data) > 1:
            corr = plot_data['Age'].corr(plot_data[col])
            correlations.append({
                'Variable': col.replace('_', ' ').title(),
                'Correlation with Age': f"{corr:.3f}",
                'Samples': len(plot_data)
            })
    
    if correlations:
        corr_df = pd.DataFrame(correlations)
        st.dataframe(corr_df, use_container_width=True)

def show_histogram_plots(df):
    """Display histogram and distribution plots."""
    
    st.subheader("📈 Distribution Analysis")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        plot_type = st.selectbox(
            "Plot type:",
            ["Age Distribution", "Telomere Length Distribution", "Mutation Rate Distribution", "Age Groups"]
        )
        
        if plot_type == "Age Groups":
            bin_size = st.slider("Age bin size (years):", 5, 20, 10)
    
    with col2:
        if plot_type == "Age Distribution":
            fig = px.histogram(
                df, 
                x='Age', 
                nbins=20,
                title="Age Distribution",
                labels={'Age': 'Age (years)', 'count': 'Number of Samples'}
            )
            st.plotly_chart(fig, use_container_width=True)
            
        elif plot_type == "Telomere Length Distribution":
            fig = px.histogram(
                df, 
                x='Telomere_Length', 
                nbins=20,
                title="Telomere Length Distribution",
                labels={'Telomere_Length': 'Telomere Length (bp)', 'count': 'Number of Samples'}
            )
            st.plotly_chart(fig, use_container_width=True)
            
        elif plot_type == "Mutation Rate Distribution":
            # Get total mutation rate column
            mutation_cols = [col for col in df.columns if 'total' in col.lower() and 'per_1k' in col]
            if mutation_cols:
                mutation_col = mutation_cols[0]
                fig = px.histogram(
                    df, 
                    x=mutation_col, 
                    nbins=20,
                    title=f"{mutation_col.replace('_', ' ').title()} Distribution",
                    labels={mutation_col: 'Mutations per 1k', 'count': 'Number of Samples'}
                )
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("No mutation rate columns found.")
                
        elif plot_type == "Age Groups":
            # Create age bins
            df_copy = df.copy()
            df_copy['Age_Group'] = pd.cut(
                df_copy['Age'], 
                bins=range(int(df['Age'].min()), int(df['Age'].max()) + bin_size, bin_size),
                include_lowest=True
            )
            
            # Count samples per age group
            age_counts = df_copy['Age_Group'].value_counts().sort_index()
            
            fig = px.bar(
                x=[str(interval) for interval in age_counts.index],
                y=age_counts.values,
                title=f"Sample Distribution by Age Groups ({bin_size}-year bins)",
                labels={'x': 'Age Group', 'y': 'Number of Samples'}
            )
            fig.update_xaxes(tickangle=45)
            st.plotly_chart(fig, use_container_width=True)

def show_custom_plots(df):
    """Display custom plot options."""
    
    st.subheader("🎨 Custom Plot Builder")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        plot_type = st.selectbox(
            "Custom plot type:",
            ["Scatter Matrix", "Correlation Heatmap", "Box Plot", "Violin Plot"]
        )
        
        if plot_type in ["Scatter Matrix", "Correlation Heatmap"]:
            # Select numeric columns
            numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
            selected_cols = st.multiselect(
                "Select columns:",
                options=numeric_cols,
                default=numeric_cols[:5] if len(numeric_cols) >= 5 else numeric_cols
            )
        else:
            # For box/violin plots
            x_var = st.selectbox("X-axis (categorical):", ['Age_Group', 'FileName'])
            y_var = st.selectbox("Y-axis (numeric):", df.select_dtypes(include=[np.number]).columns.tolist())
    
    with col2:
        if plot_type == "Scatter Matrix" and selected_cols:
            fig = px.scatter_matrix(
                df[selected_cols],
                title="Scatter Matrix of Selected Variables"
            )
            st.plotly_chart(fig, use_container_width=True)
            
        elif plot_type == "Correlation Heatmap" and selected_cols:
            corr_matrix = df[selected_cols].corr()
            fig = px.imshow(
                corr_matrix,
                text_auto=True,
                aspect="auto",
                title="Correlation Heatmap"
            )
            st.plotly_chart(fig, use_container_width=True)
            
        elif plot_type == "Box Plot":
            # Create age groups for categorical plotting
            df_copy = df.copy()
            df_copy['Age_Group'] = pd.cut(df_copy['Age'], bins=5, labels=['Young', 'Young-Adult', 'Middle', 'Adult', 'Senior'])
            
            fig = px.box(
                df_copy,
                x=x_var,
                y=y_var,
                title=f"{y_var} by {x_var}"
            )
            st.plotly_chart(fig, use_container_width=True)
            
        elif plot_type == "Violin Plot":
            # Create age groups for categorical plotting
            df_copy = df.copy()
            df_copy['Age_Group'] = pd.cut(df_copy['Age'], bins=5, labels=['Young', 'Young-Adult', 'Middle', 'Adult', 'Senior'])
            
            fig = px.violin(
                df_copy,
                x=x_var,
                y=y_var,
                title=f"{y_var} Distribution by {x_var}"
            )
            st.plotly_chart(fig, use_container_width=True)

# Call the function when the page is loaded
show_visualizations_page()
