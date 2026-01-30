import streamlit as st
import pandas as pd
import os
import tempfile
import shutil
from pathlib import Path
import sys
from datetime import datetime

# Add analysis directory to path for imports
analysis_dir = Path(__file__).parent.parent.parent / "analysis"
sys.path.append(str(analysis_dir))

# Import with explicit reload to ensure we get the latest version
import importlib
import generate_csv
importlib.reload(generate_csv)
from generate_csv import generate_csv, get_sequence_files

# Check if the function supports the callback parameter
import inspect
_generate_csv_supports_callback = 'output_callback' in inspect.signature(generate_csv).parameters
_generate_csv_supports_metadata_path = 'metadata_file_path' in inspect.signature(generate_csv).parameters

# Debug information
if st.sidebar.checkbox("Show Function Debug Info", value=False, key="function_debug"):
    st.sidebar.write(f"Generate CSV supports callback: {_generate_csv_supports_callback}")
    st.sidebar.write(f"Generate CSV supports metadata path: {_generate_csv_supports_metadata_path}")
    st.sidebar.write(f"Function signature: {inspect.signature(generate_csv)}")

def show_data_processing_page():
    """Display the data processing page for uploading and processing FASTQ/FASTA files."""
    
    st.markdown("<h2 style='text-align: center;'>📊 Data Processing</h2>", unsafe_allow_html=True)
    st.markdown("Upload and process your FASTQ/FASTA files to generate analysis-ready CSV data.")
    
    # Check for existing CSV files and show notification
    # Look in the streamlit directory (parent of pages directory)
    existing_csv_files = list(Path(__file__).parent.parent.glob("telomere_analysis*.csv"))
    
    # Debug information
    if st.sidebar.checkbox("Show CSV Debug Info", value=False, key="csv_debug"):
        st.sidebar.write(f"CSV files found: {len(existing_csv_files)}")
        st.sidebar.write(f"Data loaded: {st.session_state.get('data_loaded', False)}")
        st.sidebar.write(f"Looking in: {Path(__file__).parent.parent}")
        for f in existing_csv_files:
            st.sidebar.write(f"  - {f.name}")
        
        if st.sidebar.button("Clear Session State", help="Clear all session state data", key="clear_session"):
            for key in list(st.session_state.keys()):
                del st.session_state[key]
            st.rerun()
    
    if existing_csv_files:
        col1, col2 = st.columns([4, 1])
        with col1:
            if st.session_state.get('data_loaded', False):
                st.success(f"✅ **Found {len(existing_csv_files)} existing CSV file(s)!** Data is already loaded, but you can choose to load a different file below.")
            else:
                st.info(f"💡 **Found {len(existing_csv_files)} existing CSV file(s) in the streamlit directory!** You can choose to use these for analysis instead of processing new files.")
        with col2:
            if st.button("🔄 Refresh", help="Refresh the list of available CSV files", key="refresh_csv"):
                st.rerun()
    
    # Check if required metadata file exists
    metadata_file = Path(__file__).parent.parent.parent / "analysis" / "greider_methods_table_s2.csv"
    if not metadata_file.exists():
        st.error("⚠️ Required metadata file 'greider_methods_table_s2.csv' not found!")
        st.markdown("""
        Please ensure the metadata file is in the analysis directory. This file should contain:
        - fastq file name
        - Age (Years) 
        - Mean Telomere Length (bps)
        """)
        return
    
    # Three main options: use existing CSV, upload files, or process existing directory
    if existing_csv_files:
        processing_option = st.radio(
            "Choose processing method:",
            ["Use Existing CSV File", "Upload Files", "Process Existing Directory"],
            help="Use an existing CSV file, upload new files, or process all files in an existing directory"
        )
    else:
        processing_option = st.radio(
            "Choose processing method:",
            ["Upload Files", "Process Existing Directory"],
            help="Upload individual files or process all files in an existing directory"
        )
    
    if processing_option == "Use Existing CSV File":
        handle_existing_csv(existing_csv_files)
    elif processing_option == "Upload Files":
        handle_file_upload()
    else:
        handle_directory_processing()

def handle_existing_csv(existing_csv_files):
    """Handle using existing CSV files for analysis."""
    
    st.subheader("📄 Use Existing CSV File")
    
    if not existing_csv_files:
        st.warning("⚠️ No existing CSV files found in the streamlit directory.")
        return
    
    # Display available CSV files
    st.info(f"Found {len(existing_csv_files)} existing CSV file(s):")
    
    # Show file details
    csv_data = []
    for i, csv_file in enumerate(existing_csv_files):
        try:
            # Get file stats
            file_size = csv_file.stat().st_size
            file_mtime = datetime.fromtimestamp(csv_file.stat().st_mtime)
            
            # Load a preview of the data
            df_preview = pd.read_csv(csv_file, nrows=5)
            total_rows = len(pd.read_csv(csv_file))
            
            csv_data.append({
                'file': csv_file,
                'name': csv_file.name,
                'size': file_size,
                'modified': file_mtime,
                'rows': total_rows,
                'columns': len(df_preview.columns),
                'preview': df_preview
            })
            
        except Exception as e:
            st.error(f"Error reading {csv_file.name}: {e}")
            continue
    
    # Display file selection
    if csv_data:
        # Create a selection interface
        selected_file = None
        
        if len(csv_data) == 1:
            selected_file = csv_data[0]
            st.success(f"✅ Using the only available file: {selected_file['name']}")
        else:
            # Multiple files - let user choose
            file_options = [f"{f['name']} ({f['rows']} rows, {f['size']:,} bytes, modified: {f['modified'].strftime('%Y-%m-%d %H:%M')})" 
                          for f in csv_data]
            selected_index = st.selectbox(
                "Select a CSV file to use:",
                range(len(csv_data)),
                format_func=lambda x: file_options[x],
                help="Choose which CSV file to load for analysis"
            )
            selected_file = csv_data[selected_index]
        
        if selected_file:
            # Show file information
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("File Size", f"{selected_file['size']:,} bytes")
            
            with col2:
                st.metric("Total Rows", selected_file['rows'])
            
            with col3:
                st.metric("Total Columns", selected_file['columns'])
            
            with col4:
                st.metric("Last Modified", selected_file['modified'].strftime('%Y-%m-%d'))
            
            # Show data preview
            with st.expander("📋 Data Preview", expanded=True):
                st.dataframe(selected_file['preview'], use_container_width=True)
            
            # Show column information
            with st.expander("📊 Column Information"):
                try:
                    df_full = pd.read_csv(selected_file['file'])
                    col_info = pd.DataFrame({
                        'Column': df_full.columns,
                        'Type': df_full.dtypes,
                        'Non-Null Count': df_full.count(),
                        'Null Count': df_full.isnull().sum()
                    })
                    st.dataframe(col_info, use_container_width=True)
                except Exception as e:
                    st.error(f"Error loading full data: {e}")
            
            # Load data button
            if st.button("🔄 Load This CSV File", type="primary", key="load_csv"):
                try:
                    # Load the full dataset
                    df = pd.read_csv(selected_file['file'])
                    
                    # Update session state
                    st.session_state.csv_generated = True
                    st.session_state.data_loaded = True
                    st.session_state.analysis_data = df
                    
                    st.success("✅ CSV file loaded successfully!")
                    st.info("💡 Data is now ready for analysis! Navigate to other pages to explore your data.")
                    
                    # Show summary statistics
                    show_data_summary(df)
                    
                except Exception as e:
                    st.error(f"❌ Error loading CSV file: {e}")

def show_data_summary(df):
    """Display summary statistics for the loaded data."""
    
    st.subheader("📊 Data Summary")
    
    # Basic statistics
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Total Samples", len(df))
    
    with col2:
        st.metric("Total Features", len(df.columns))
    
    with col3:
        if 'Age' in df.columns:
            age_range = f"{df['Age'].min():.0f} - {df['Age'].max():.0f}"
            st.metric("Age Range (years)", age_range)
        else:
            st.metric("Age Column", "Not found")
    
    # Data quality metrics
    st.subheader("📈 Data Quality")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        missing_data = df.isnull().sum().sum()
        st.metric("Missing Values", missing_data)
    
    with col2:
        duplicate_rows = df.duplicated().sum()
        st.metric("Duplicate Rows", duplicate_rows)
    
    with col3:
        if 'Total_Reads' in df.columns:
            avg_reads = df['Total_Reads'].mean()
            st.metric("Avg Total Reads", f"{avg_reads:,.0f}")
        else:
            st.metric("Total Reads", "Not found")
    
    with col4:
        if 'Telomere_Length' in df.columns:
            avg_length = df['Telomere_Length'].mean()
            st.metric("Avg Telomere Length", f"{avg_length:,.0f} bps")
        else:
            st.metric("Telomere Length", "Not found")

def handle_file_upload():
    """Handle individual file uploads."""
    
    st.subheader("📁 Upload FASTQ/FASTA Files")
    
    uploaded_files = st.file_uploader(
        "Choose FASTQ or FASTA files",
        type=['fastq', 'fasta', 'fa', 'fas', 'gz'],
        accept_multiple_files=True,
        help="Upload one or more sequencing files. Supports .fastq, .fasta, .fa, .fas, and compressed .gz files"
    )
    
    if uploaded_files:
        st.success(f"✅ {len(uploaded_files)} files uploaded successfully!")
        
        # Show uploaded files
        with st.expander("📋 View Uploaded Files"):
            for i, file in enumerate(uploaded_files):
                st.write(f"{i+1}. {file.name} ({file.size:,} bytes)")
        
        # Process files
        if st.button("🔄 Process Uploaded Files", type="primary", key="process_uploaded"):
            # Add console output section
            st.subheader("🖥️ Processing Console Output")
            st.info("Real-time processing output will appear below:")
            
            # Create a container for console output
            with st.container():
                process_uploaded_files(uploaded_files)

def handle_directory_processing():
    """Handle processing files from an existing directory."""
    
    st.subheader("📂 Process Existing Directory")
    
    # Directory input
    data_dir = st.text_input(
        "Enter directory path containing FASTQ/FASTA files:",
        value="/Users/akhilpeddikuppa/FieldLab/GreiderCodeSearch/greider_data_download",
        help="Path to directory containing your sequencing files"
    )
    
    if data_dir and os.path.exists(data_dir):
        # Check for sequence files
        try:
            sequence_files = get_sequence_files(data_dir)
            if sequence_files:
                st.success(f"✅ Found {len(sequence_files)} sequence files in directory")
                
                # Show found files
                with st.expander("📋 View Found Files"):
                    for i, file_path in enumerate(sequence_files[:10]):  # Show first 10
                        st.write(f"{i+1}. {os.path.basename(file_path)}")
                    if len(sequence_files) > 10:
                        st.write(f"... and {len(sequence_files) - 10} more files")
                
                # Process files
                if st.button("🔄 Process Directory Files", type="primary", key="process_directory"):
                    # Add console output section
                    st.subheader("🖥️ Processing Console Output")
                    st.info("Real-time processing output will appear below:")
                    
                    # Create a container for console output
                    with st.container():
                        process_directory_files(data_dir)
            else:
                st.warning("⚠️ No FASTQ/FASTA files found in the specified directory")
        except Exception as e:
            st.error(f"Error scanning directory: {e}")
    elif data_dir:
        st.error("❌ Directory does not exist. Please check the path.")

def process_uploaded_files(uploaded_files):
    """Process uploaded files by saving them temporarily and running analysis."""
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    # Create console output display
    console_output = st.empty()
    console_messages = []
    
    def update_console(message):
        """Update console output display"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        formatted_message = f"[{timestamp}] {message}"
        console_messages.append(formatted_message)
        
        # Create a scrollable text area with better formatting
        console_text = '\n'.join(console_messages[-30:])  # Show last 30 lines
        console_output.text_area(
            "Console Output", 
            value=console_text, 
            height=400, 
            disabled=True,
            help="Real-time processing output. Scroll to see all messages."
        )
    
    try:
        # Create temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            status_text.text("💾 Saving uploaded files...")
            progress_bar.progress(0.1)
            update_console("💾 Saving uploaded files...")
            
            # Save uploaded files to temporary directory
            temp_files = []
            for i, uploaded_file in enumerate(uploaded_files):
                temp_path = os.path.join(temp_dir, uploaded_file.name)
                with open(temp_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())
                temp_files.append(temp_path)
                progress_bar.progress(0.1 + (i + 1) / len(uploaded_files) * 0.3)
                update_console(f"Saved: {uploaded_file.name}")
            
            status_text.text("🔄 Processing files...")
            progress_bar.progress(0.4)
            update_console("🔄 Starting file processing...")
            
            # Generate CSV using the existing function with console output callback
            metadata_file = Path(__file__).parent.parent.parent / "analysis" / "greider_methods_table_s2.csv"
            
            # Build function call arguments based on what the function supports
            kwargs = {}
            if _generate_csv_supports_callback:
                kwargs['output_callback'] = update_console
            if _generate_csv_supports_metadata_path:
                kwargs['metadata_file_path'] = str(metadata_file)
            
            if not _generate_csv_supports_callback:
                update_console("⚠️ Using fallback mode (console output limited)")
            
            generate_csv(temp_dir, **kwargs)
            
            progress_bar.progress(1.0)
            status_text.text("✅ Processing complete!")
            update_console("✅ Processing complete!")
            
            # Update session state
            st.session_state.csv_generated = True
            
            # Show results
            show_processing_results()
            
    except Exception as e:
        st.error(f"❌ Error processing files: {e}")
        progress_bar.progress(0)
        status_text.text("❌ Processing failed")
        update_console(f"❌ Error: {e}")

def process_directory_files(data_dir):
    """Process files from an existing directory."""
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    # Create console output display
    console_output = st.empty()
    console_messages = []
    
    def update_console(message):
        """Update console output display"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        formatted_message = f"[{timestamp}] {message}"
        console_messages.append(formatted_message)
        
        # Create a scrollable text area with better formatting
        console_text = '\n'.join(console_messages[-30:])  # Show last 30 lines
        console_output.text_area(
            "Console Output", 
            value=console_text, 
            height=400, 
            disabled=True,
            help="Real-time processing output. Scroll to see all messages."
        )
    
    try:
        status_text.text("🔄 Processing files from directory...")
        progress_bar.progress(0.5)
        update_console(f"🔄 Processing files from directory: {data_dir}")
        
        # Generate CSV using the existing function with console output callback
        metadata_file = Path(__file__).parent.parent.parent / "analysis" / "greider_methods_table_s2.csv"
        
        # Build function call arguments based on what the function supports
        kwargs = {}
        if _generate_csv_supports_callback:
            kwargs['output_callback'] = update_console
        if _generate_csv_supports_metadata_path:
            kwargs['metadata_file_path'] = str(metadata_file)
        
        if not _generate_csv_supports_callback:
            update_console("⚠️ Using fallback mode (console output limited)")
        
        generate_csv(data_dir, **kwargs)
        
        progress_bar.progress(1.0)
        status_text.text("✅ Processing complete!")
        update_console("✅ Processing complete!")
        
        # Update session state
        st.session_state.csv_generated = True
        
        # Show results
        show_processing_results()
        
    except Exception as e:
        st.error(f"❌ Error processing files: {e}")
        progress_bar.progress(0)
        status_text.text("❌ Processing failed")
        update_console(f"❌ Error: {e}")

def show_processing_results():
    """Display the results of CSV generation."""
    
    st.subheader("📊 Processing Results")
    
    # Check for generated CSV files
    csv_files = list(Path('.').glob('telomere_analysis*.csv'))
    
    if csv_files:
        st.success("✅ CSV file generated successfully!")
        
        # Load and display data preview
        try:
            df = pd.read_csv(csv_files[0])
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.metric("Total Samples", len(df))
            
            with col2:
                st.metric("Total Features", len(df.columns))
            
            with col3:
                age_range = f"{df['Age'].min():.0f} - {df['Age'].max():.0f}"
                st.metric("Age Range (years)", age_range)
            
            # Data preview
            st.subheader("📋 Data Preview")
            st.dataframe(df.head(10), use_container_width=True)
            
            # Column information
            with st.expander("📊 Column Information"):
                col_info = pd.DataFrame({
                    'Column': df.columns,
                    'Type': df.dtypes,
                    'Non-Null Count': df.count(),
                    'Null Count': df.isnull().sum()
                })
                st.dataframe(col_info, use_container_width=True)
            
            # Update session state
            st.session_state.data_loaded = True
            st.session_state.analysis_data = df
            
            st.info("💡 Data is now ready for analysis! Navigate to other pages to explore your data.")
            
            # Show summary statistics
            show_data_summary(df)
            
        except Exception as e:
            st.error(f"Error loading generated CSV: {e}")
    else:
        st.error("❌ No CSV file was generated. Please check the error messages above.")

def show_processing_info():
    """Display information about the processing pipeline."""
    
    st.subheader("ℹ️ Processing Information")
    
    st.markdown("""
    ### What happens during processing:
    
    1. **File Reading**: FASTQ/FASTA files are read and sequences extracted
    2. **Pattern Matching**: Telomere sequences and mutations are identified using predefined patterns
    3. **Counting**: Raw counts are calculated for each pattern type
    4. **Normalization**: Counts are normalized per 1000 bases for comparison
    5. **Feature Engineering**: Additional features like composite scores are calculated
    6. **Metadata Integration**: Age and telomere length data are merged from metadata file
    
    ### Supported File Formats:
    - `.fastq` - Standard FASTQ format
    - `.fasta` / `.fa` / `.fas` - FASTA format
    - `.gz` - Compressed versions of above formats
    
    ### Required Metadata:
    The `greider_methods_table_s2.csv` file must contain:
    - `fastq file name` - Matching filename (without extension)
    - `Age (Years)` - Age of the sample
    - `Mean Telomere Length (bps)` - Telomere length measurement
    """)
    
    # Show current metadata file info
    metadata_file = Path(__file__).parent.parent.parent / "analysis" / "greider_methods_table_s2.csv"
    if metadata_file.exists():
        try:
            metadata_df = pd.read_csv(metadata_file)
            st.success(f"✅ Metadata file loaded: {len(metadata_df)} samples")
            
            with st.expander("📋 Metadata Preview"):
                st.dataframe(metadata_df.head(), use_container_width=True)
        except Exception as e:
            st.error(f"Error reading metadata file: {e}")
    else:
        st.warning("⚠️ Metadata file not found")

# Call the function when the page is loaded
show_data_processing_page()
