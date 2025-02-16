import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Create directories if they don't exist
for dir_path in ['proteomics_analysis/data', 'proteomics_analysis/plots', 'proteomics_analysis/results']:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

print("Loading and preprocessing data files...")

# Read the Excel file
excel_file = "1-s2.0-S1097276521009540-mmc6.xlsx"

# Read all sheets with correct names
update_df = pd.read_excel(excel_file, sheet_name='Uptake')
translation_df = pd.read_excel(excel_file, sheet_name='Translation')

# Convert binary columns to numeric
binary_columns = [
    'Significantly decreased (CCCP/DMSO)',
    'Uptake defect',
    'Mitochondrial small ribosomal subunit',
    'Mitochondrial genome-encoded proteins',
    'Respiratory complex I',
    'TOM complex',
    'TIM22 TIM23 PAM complexes'
]

# Function to convert values to binary
def to_binary(val):
    if pd.isna(val):
        return 0
    if isinstance(val, (int, float)):
        return 1 if val >= 1 else 0
    if isinstance(val, str):
        return 1 if val.lower() in ['yes', 'true', '1'] else 0
    return 0

# Convert binary columns
for col in binary_columns:
    if col in update_df.columns:
        update_df[col] = update_df[col].apply(to_binary)

print("\nColumn value counts for binary columns:")
for col in binary_columns:
    if col in update_df.columns:
        print(f"\n{col}:")
        print(update_df[col].value_counts())

def analyze_distribution(df, column_name, save_path):
    """Analyze and plot distribution of a column"""
    plt.figure(figsize=(10, 6))
    sns.histplot(data=df, x=column_name, kde=True)
    plt.title(f'Distribution of {column_name}')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_volcano_plot(df, fc_column, pval_column, save_path, title):
    """Create volcano plot"""
    plt.figure(figsize=(12, 8))
    
    # Create scatter plot
    significant = df[pval_column] < 0.05
    
    # Plot non-significant points
    plt.scatter(df[~significant][fc_column], 
               -np.log10(df[~significant][pval_column]), 
               alpha=0.5, 
               color='gray',
               label='Not Significant')
    
    # Plot significant points
    plt.scatter(df[significant][fc_column], 
               -np.log10(df[significant][pval_column]), 
               alpha=0.7,
               color='red',
               label='Significant (p<0.05)')
    
    # Add threshold lines
    plt.axhline(y=-np.log10(0.05), color='r', linestyle='--', alpha=0.3)
    plt.axvline(x=-1, color='r', linestyle='--', alpha=0.3)
    plt.axvline(x=1, color='r', linestyle='--', alpha=0.3)
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10(Adjusted P-value)')
    plt.title(title)
    plt.legend()
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

print("\nAnalyzing Uptake sheet (Mitochondrial Import)...")

print("\nCreating volcano plots...")
# Create volcano plots for CCCP vs DMSO
create_volcano_plot(
    update_df,
    'Log2(CCCP/DMSO)',
    'Adjusted P value (q value; CCCP vs DMSO)',
    'proteomics_analysis/plots/volcano_cccp_vs_dmso_update.png',
    'Volcano Plot: CCCP vs DMSO (Uptake)'
)

# Create volcano plots for CCCP+ISRIB vs DMSO
create_volcano_plot(
    update_df,
    'Log2(CCCP+ISRIB/DMSO)',
    'Adjusted P value (q value; CCCP+ISRIB vs DMSO)',
    'proteomics_analysis/plots/volcano_cccp_isrib_vs_dmso_update.png',
    'Volcano Plot: CCCP+ISRIB vs DMSO (Uptake)'
)

# Distribution analysis
analyze_distribution(
    update_df,
    'Log2(CCCP/DMSO)',
    'proteomics_analysis/plots/distribution_cccp_dmso_update.png'
)

print("\nAnalyzing protein categories...")
# Create summary statistics for different protein categories
protein_categories = {
    'Mitochondrial small ribosomal subunit': 'Mito Ribosome',
    'Mitochondrial genome-encoded proteins': 'Mt-encoded',
    'Respiratory complex I': 'Complex I',
    'TOM complex': 'TOM',
    'TIM22 TIM23 PAM complexes': 'TIM/PAM'
}

# Save basic statistics and category analysis
with open('proteomics_analysis/results/initial_analysis.txt', 'w') as f:
    f.write("Initial Analysis Results\n")
    f.write("=======================\n\n")
    
    # Uptake sheet statistics
    f.write("Uptake Sheet Statistics:\n")
    f.write(f"Total proteins analyzed: {len(update_df)}\n")
    
    # Count proteins with binary indicators
    sig_decreased = update_df['Significantly decreased (CCCP/DMSO)'].sum()
    uptake_defect = update_df['Uptake defect'].sum()
    
    f.write(f"Significantly decreased proteins (CCCP/DMSO): {int(sig_decreased)}\n")
    f.write(f"Proteins with uptake defect: {int(uptake_defect)}\n")
    
    # Add more detailed statistics
    f.write("\nStatistical Summary:\n")
    f.write("\nCCCP vs DMSO Summary:\n")
    f.write(f"Mean Log2FC: {update_df['Log2(CCCP/DMSO)'].mean():.3f}\n")
    f.write(f"Median Log2FC: {update_df['Log2(CCCP/DMSO)'].median():.3f}\n")
    f.write(f"Std Dev Log2FC: {update_df['Log2(CCCP/DMSO)'].std():.3f}\n")
    
    significant_proteins = update_df[update_df['Adjusted P value (q value; CCCP vs DMSO)'] < 0.05]
    f.write(f"\nTotal significant proteins (p<0.05): {len(significant_proteins)}\n")
    f.write(f"Significant upregulated (Log2FC > 1): {len(significant_proteins[significant_proteins['Log2(CCCP/DMSO)'] > 1])}\n")
    f.write(f"Significant downregulated (Log2FC < -1): {len(significant_proteins[significant_proteins['Log2(CCCP/DMSO)'] < -1])}\n")
    
    # Category analysis
    f.write("\nProtein Category Analysis:\n")
    for col, name in protein_categories.items():
        if col in update_df.columns:
            category_proteins = update_df[update_df[col] == 1]
            sig_category = category_proteins[category_proteins['Adjusted P value (q value; CCCP vs DMSO)'] < 0.05]
            
            if len(category_proteins) > 0:
                f.write(f"\n{name}:\n")
                f.write(f"Total proteins: {len(category_proteins)}\n")
                f.write(f"Significant proteins: {len(sig_category)}\n")
                f.write(f"Mean Log2FC: {category_proteins['Log2(CCCP/DMSO)'].mean():.3f}\n")
                
                # Additional category-specific analysis
                decreased = len(category_proteins[category_proteins['Log2(CCCP/DMSO)'] < -1])
                increased = len(category_proteins[category_proteins['Log2(CCCP/DMSO)'] > 1])
                f.write(f"Decreased (Log2FC < -1): {decreased}\n")
                f.write(f"Increased (Log2FC > 1): {increased}\n")

print("\nInitial analysis complete. Check the 'proteomics_analysis' directory for results.")

# Additional analysis of translation data
print("\nAnalyzing Translation sheet...")
if not translation_df.empty:
    create_volcano_plot(
        translation_df,
        'Log2(CCCP/DMSO)',
        'Adjusted P value (q value; CCCP vs DMSO)',
        'proteomics_analysis/plots/volcano_cccp_vs_dmso_translation.png',
        'Volcano Plot: CCCP vs DMSO (Translation)'
    )
