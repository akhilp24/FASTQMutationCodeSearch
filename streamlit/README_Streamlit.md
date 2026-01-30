# Telomere Analysis Dashboard

An interactive Streamlit web application for comprehensive telomere sequencing data analysis.

## Features

### 📊 Data Processing

- Upload FASTQ/FASTA files or process existing directories
- Generate analysis-ready CSV files from sequencing data
- Real-time progress tracking and data validation
- Integration with metadata files for age and telomere length

### 📈 Interactive Visualizations

- **Scatter Plots**: Interactive plots with customizable variables, colors, and trendlines
- **Trendline Analysis**: Correlation analysis with R² calculations
- **Distribution Plots**: Histograms, box plots, and violin plots
- **Custom Plots**: Scatter matrices, correlation heatmaps, and more

### 📊 Statistical Analysis

- **Spearman Correlations**: Comprehensive correlation analysis with age
- **Composite Scores**: Advanced feature engineering and scoring
- **Correlation Heatmaps**: Visualize pairwise relationships
- **Curve Fitting**: Multiple curve types (linear, exponential, logarithmic, polynomial, power)

### 🧬 Mutational Signatures

- **COSMIC SBS Signatures**: SBS1, SBS4, SBS5, and SBS18 analysis
- **Signature Correlations**: Age correlation analysis for each signature
- **Interactive Plots**: Detailed signature visualization with trendlines

### 💾 Export & Reporting

- **Data Export**: CSV downloads for filtered and full datasets
- **Plot Export**: High-resolution PNG, PDF, and SVG downloads
- **Report Generation**: Comprehensive HTML, PDF, and Markdown reports
- **Batch Download**: Single ZIP file with all results

## Installation

1. **Install Dependencies**:

   ```bash
   pip install -r requirements.txt
   ```

2. **Required Files**:
   - Ensure `greider_methods_table_s2.csv` is in the analysis directory
   - This file should contain:
     - `fastq file name`: Matching filename (without extension)
     - `Age (Years)`: Age of the sample
     - `Mean Telomere Length (bps)`: Telomere length measurement

## Usage

1. **Start the Application**:

   ```bash
   streamlit run app.py
   ```

2. **Navigate the Dashboard**:
   - **Home**: Project overview and quick stats
   - **Data Processing**: Upload and process your sequencing files
   - **Visualizations**: Explore data through interactive plots
   - **Statistical Analysis**: Perform comprehensive statistical analysis
   - **Mutational Signatures**: Analyze COSMIC SBS signatures
   - **Export Results**: Download plots, data, and reports

## File Structure

```
analysis/
├── app.py                          # Main Streamlit application
├── pages/                          # Page modules
│   ├── data_processing.py          # File upload and CSV generation
│   ├── visualizations.py           # Interactive plotting
│   ├── statistical_analysis.py     # Statistical analysis tools
│   ├── mutational_signatures.py    # SBS signature analysis
│   └── export_results.py           # Export and download tools
├── requirements.txt                # Python dependencies
├── README_Streamlit.md             # This file
├── generate_csv.py                 # Original CSV generation module
├── plotting.py                     # Original plotting functions
├── scatterplot.py                  # Original scatter plot functions
├── trendline.py                    # Original trendline functions
├── histogram.py                    # Original histogram functions
├── signatures.py                   # Original signature analysis
└── main.py                         # Original main script
```

## Key Features

### Interactive Controls

- **Age Range Sliders**: Filter data by age ranges
- **Sample Selection**: Choose specific samples for analysis
- **Variable Selection**: Customize plots with different variables
- **Real-time Updates**: Dynamic plot updates based on selections

### Data Processing Pipeline

1. **File Reading**: Supports FASTQ, FASTA, and compressed files
2. **Pattern Matching**: Identifies telomere sequences and mutations
3. **Normalization**: Per-1000-base normalization for comparison
4. **Feature Engineering**: Composite scores and advanced metrics
5. **Metadata Integration**: Age and telomere length from reference file

### Analysis Methods

- **Strand-specific Analysis**: Separate C-strand and G-strand analysis
- **Mutation Counting**: Comprehensive mutation pattern detection
- **Statistical Testing**: Spearman correlations and significance testing
- **Curve Fitting**: Multiple mathematical models for relationship analysis
- **Signature Analysis**: COSMIC mutational signature scoring

## Supported File Formats

- **Input**: `.fastq`, `.fasta`, `.fa`, `.fas`, `.gz` (compressed)
- **Output**: `.csv`, `.png`, `.pdf`, `.svg`, `.html`, `.md`, `.zip`

## Browser Compatibility

- Chrome (recommended)
- Firefox
- Safari
- Edge

## Performance Notes

- **Large Files**: Processing large FASTQ files may take several minutes
- **Memory Usage**: Large datasets may require significant RAM
- **Caching**: Results are cached in session state for faster navigation
- **Progress Tracking**: Real-time progress bars for long operations

## Troubleshooting

### Common Issues

1. **"No data loaded" Warning**:

   - Go to Data Processing page first
   - Upload files or process existing directory
   - Ensure metadata file is present

2. **File Upload Errors**:

   - Check file format (FASTQ/FASTA)
   - Ensure files are not corrupted
   - Verify file permissions

3. **Missing Metadata**:

   - Ensure `greider_methods_table_s2.csv` exists
   - Check filename matching between data and metadata
   - Verify CSV format and column names

4. **Plot Display Issues**:
   - Refresh the page
   - Check browser compatibility
   - Clear browser cache

### Performance Optimization

- **Filter Data**: Use age range and sample selection to reduce data size
- **Close Tabs**: Close unused browser tabs to free memory
- **Restart App**: Restart Streamlit if memory usage is high

## Contributing

To extend the dashboard:

1. **Add New Pages**: Create new modules in the `pages/` directory
2. **Modify Analysis**: Update existing analysis functions
3. **Add Visualizations**: Extend plotting capabilities
4. **Improve UI**: Enhance user interface and experience

## License

This project is part of the Greider Telomere Analysis research project.

## Support

For technical support or questions about the analysis methods, please refer to the original analysis scripts or contact the research team.
