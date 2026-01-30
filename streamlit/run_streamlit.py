#!/usr/bin/env python3
"""
Telomere Analysis Dashboard - Streamlit Application Launcher

Usage:
    python run_streamlit.py

Or directly:
    streamlit run app.py
"""

import subprocess
import sys
import os
from pathlib import Path

def check_dependencies():
    """Check if required packages are installed."""
    required_packages = [
        'streamlit',
        'pandas', 
        'numpy',
        'matplotlib',
        'seaborn',
        'plotly',
        'scipy'
    ]
    
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        print("❌ Missing required packages:")
        for package in missing_packages:
            print(f"   - {package}")
        print("\n📦 Install missing packages with:")
        print("   pip install -r requirements.txt")
        return False
    
    print("✅ All required packages are installed!")
    return True

def check_metadata_file():
    """Check if required metadata file exists."""
    analysis_dir = Path(__file__).parent.parent / "analysis"
    metadata_file = analysis_dir / "greider_methods_table_s2.csv"
    
    if not metadata_file.exists():
        print("⚠️  Warning: Required metadata file 'greider_methods_table_s2.csv' not found!")
        print("   This file should contain:")
        print("   - fastq file name")
        print("   - Age (Years)")
        print("   - Mean Telomere Length (bps)")
        print("   Please ensure this file is in the analysis directory.")
        return False
    
    print("✅ Metadata file found!")
    return True

def main():
    """Main launcher function."""
    print("🧬 Telomere Analysis Dashboard")
    print("=" * 40)
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Check metadata file
    check_metadata_file()
    
    print("\n🚀 Starting Streamlit application...")
    print("   The dashboard will open in your default web browser.")
    print("   Press Ctrl+C to stop the application.")
    print("=" * 40)
    
    try:
        # Launch Streamlit
        subprocess.run([
            sys.executable, "-m", "streamlit", "run", "app.py",
            "--server.port", "8501",
            "--server.address", "localhost"
        ])
    except KeyboardInterrupt:
        print("\n👋 Application stopped by user.")
    except Exception as e:
        print(f"\n❌ Error launching application: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
