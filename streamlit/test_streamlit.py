#!/usr/bin/env python3
"""
Test script for the Telomere Analysis Dashboard

This script performs basic tests to ensure the Streamlit app components work correctly.
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np

# Add current directory to path
sys.path.append(str(Path(__file__).parent))

def test_imports():
    """Test if all required modules can be imported."""
    print("🧪 Testing imports...")
    
    try:
        import streamlit as st
        print("✅ Streamlit imported successfully")
    except ImportError as e:
        print(f"❌ Streamlit import failed: {e}")
        return False
    
    try:
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns
        import plotly.express as px
        import plotly.graph_objects as go
        from scipy import stats
        print("✅ All analysis libraries imported successfully")
    except ImportError as e:
        print(f"❌ Analysis library import failed: {e}")
        return False
    
    try:
        from pages import data_processing, visualizations, statistical_analysis, mutational_signatures, export_results
        print("✅ All page modules imported successfully")
    except ImportError as e:
        print(f"❌ Page module import failed: {e}")
        return False
    
    return True

def test_data_processing():
    """Test data processing functions."""
    print("\n🧪 Testing data processing...")
    
    try:
        # Add analysis directory to path
        analysis_dir = Path(__file__).parent.parent / "analysis"
        sys.path.append(str(analysis_dir))
        from generate_csv import get_sequence_files
        print("✅ CSV generation functions imported")
    except ImportError as e:
        print(f"❌ CSV generation import failed: {e}")
        return False
    
    return True

def test_sample_data():
    """Test with sample data."""
    print("\n🧪 Testing with sample data...")
    
    # Create sample data
    sample_data = pd.DataFrame({
        'FileName': ['sample1.fastq', 'sample2.fastq', 'sample3.fastq'],
        'Age': [25, 45, 65],
        'Telomere_Length': [8000, 6000, 4000],
        'Total_Reads': [10000, 12000, 8000],
        'g_strand_mutations_T>C_t1_per_1k': [0.5, 1.2, 2.1],
        'c_strand_mutations_C>T_c1_per_1k': [0.3, 0.8, 1.5]
    })
    
    print(f"✅ Sample data created with {len(sample_data)} samples")
    
    # Test basic operations
    try:
        # Test correlation
        from scipy import stats
        corr, p_val = stats.spearmanr(sample_data['Age'], sample_data['g_strand_mutations_T>C_t1_per_1k'])
        print(f"✅ Correlation test passed: r={corr:.3f}, p={p_val:.3f}")
        
        # Test plotting
        import plotly.express as px
        fig = px.scatter(sample_data, x='Age', y='g_strand_mutations_T>C_t1_per_1k')
        print("✅ Plotly plotting test passed")
        
    except Exception as e:
        print(f"❌ Sample data test failed: {e}")
        return False
    
    return True

def test_file_structure():
    """Test if all required files exist."""
    print("\n🧪 Testing file structure...")
    
    required_files = [
        'app.py',
        'pages/data_processing.py',
        'pages/visualizations.py', 
        'pages/statistical_analysis.py',
        'pages/mutational_signatures.py',
        'pages/export_results.py',
        'requirements.txt',
        'README_Streamlit.md'
    ]
    
    missing_files = []
    for file in required_files:
        if not Path(file).exists():
            missing_files.append(file)
    
    if missing_files:
        print("❌ Missing files:")
        for file in missing_files:
            print(f"   - {file}")
        return False
    
    print("✅ All required files present")
    return True

def main():
    """Run all tests."""
    print("🧬 Telomere Analysis Dashboard - Test Suite")
    print("=" * 50)
    
    tests = [
        test_imports,
        test_data_processing,
        test_sample_data,
        test_file_structure
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
        else:
            print("❌ Test failed!")
    
    print("\n" + "=" * 50)
    print(f"📊 Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! The Streamlit app should work correctly.")
        print("\n🚀 To start the app, run:")
        print("   python run_streamlit.py")
        print("   or")
        print("   streamlit run app.py")
    else:
        print("⚠️  Some tests failed. Please check the errors above.")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
