# Mitochondrial Proteomics Analysis

## Overview

This analysis examines TMT-based mePROD proteomics data investigating mitochondrial protein import and translation under CCCP treatment and ISRIB rescue conditions.

## Directory Structure

```
proteomics_analysis/
├── data/            # Data directory
├── plots/           # Generated visualizations
│   ├── category_rescue_percentages.png
│   ├── category_treatment_effects.png
│   ├── category_rescue_magnitudes.png
│   ├── category_metrics_heatmap.png
│   ├── volcano_cccp_vs_dmso_update.png
│   └── volcano_cccp_isrib_vs_dmso_update.png
├── results/         # Analysis results
│   ├── initial_analysis.txt         # Basic statistical analysis
│   ├── biological_insights.txt      # Primary biological findings
│   ├── isrib_rescue_insights.txt    # ISRIB-specific analysis
│   ├── category_rescue_insights.txt # Category-specific patterns
│   ├── category_analysis_summary.txt # Detailed category statistics
│   └── final_summary.txt           # Comprehensive analysis summary
└── *.py            # Analysis scripts
```

## Key Files

### Analysis Scripts

- `proteomics_analysis.py`: Initial data processing and analysis
- `detailed_analysis.py`: In-depth analysis of rescue effects
- `category_visualization.py`: Category-specific visualizations

### Results

1. **Initial Analysis** (`initial_analysis.txt`)

   - Basic statistics
   - Protein counts
   - Global effects

2. **Biological Insights** (`biological_insights.txt`)

   - Primary findings
   - Initial hypotheses
   - Pattern analysis

3. **ISRIB Rescue Analysis** (`isrib_rescue_insights.txt`)

   - ISRIB-specific effects
   - Rescue patterns
   - Mechanism insights

4. **Category Analysis** (`category_rescue_insights.txt`)

   - Category-specific patterns
   - Comparative analysis
   - Hierarchical effects

5. **Final Summary** (`final_summary.txt`)
   - Comprehensive overview
   - Key findings
   - Future directions

## Key Findings

### Global Effects

- Total proteins analyzed: 4608
- CCCP causes widespread protein decrease
- ISRIB rescues 2423 proteins
- Mean rescue magnitude: 1.536 Log2FC

### Category-Specific Patterns

1. **TOM Complex** (Import Machinery)

   - 83.3% rescue rate
   - Near-baseline recovery
   - Strongest rescue magnitude

2. **Complex I** (Respiratory Chain)

   - 67.3% rescue rate
   - Substantial recovery
   - Coordinated rescue pattern

3. **Mitochondrial Ribosomes**

   - 61.5% rescue rate
   - Moderate recovery
   - Intermediate rescue effects

4. **Mt-encoded Proteins**
   - 28.6% rescue rate
   - Negative rescue effect
   - Unique vulnerability pattern

## Visualizations

- Category rescue percentages
- Treatment effects comparison
- Rescue magnitude analysis
- Normalized metrics heatmap
- Volcano plots for treatment comparisons

## Usage

To rerun the analysis:

1. Ensure required Python packages are installed
2. Run scripts in order:
   ```bash
   python proteomics_analysis.py
   python detailed_analysis.py
   python category_visualization.py
   ```
3. Check the plots/ and results/ directories for outputs

## Data Source

TMT-based mePROD proteomics data examining:

- CCCP vs DMSO comparison
- CCCP+ISRIB vs DMSO comparison
- Focus on mitochondrial proteins
- Multiple protein categories

## Contact

For questions or feedback about this analysis, please refer to the original publication or contact the authors.
