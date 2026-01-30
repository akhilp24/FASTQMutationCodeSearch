"""
Utility functions for the Streamlit telomere analysis app.
"""

import streamlit as st
import pandas as pd
from pathlib import Path

def load_data_if_available():
    """
    Load data from session state or automatically load existing CSV file.
    
    Returns:
        tuple: (df, success_message) where df is the loaded DataFrame and 
               success_message is a string to display to the user
    """
    # Check if data is already loaded in session state
    if st.session_state.get('data_loaded', False):
        return st.session_state.analysis_data, None
    
    # Try to load existing CSV file automatically
    existing_csv_files = list(Path(__file__).parent.glob("telomere_analysis*.csv"))
    
    if existing_csv_files:
        try:
            df = pd.read_csv(existing_csv_files[0])
            st.session_state.data_loaded = True
            st.session_state.analysis_data = df
            success_msg = f"✅ Loaded {existing_csv_files[0].name} with {len(df)} samples and {len(df.columns)} features!"
            return df, success_msg
        except Exception as e:
            error_msg = f"❌ Error loading CSV file: {e}"
            return None, error_msg
    else:
        return None, "⚠️ No data loaded. Please go to the Data Processing page first."

def check_data_loaded():
    """
    Check if data is loaded and display appropriate messages.
    
    Returns:
        bool: True if data is loaded, False otherwise
    """
    df, message = load_data_if_available()
    
    if df is not None:
        if message:
            st.success(message)
        return True
    else:
        if "Error loading" in message:
            st.error(message)
            st.warning("⚠️ Please go to the Data Processing page to load data manually.")
        else:
            st.warning(message)
        return False

