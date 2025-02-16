import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read the data
df = pd.read_excel("1-s2.0-S1097276521009540-mmc6.xlsx", sheet_name='Uptake')

# Convert binary columns more carefully
binary_columns = [
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
        val = val.lower().strip()
        return 1 if val in ['yes', 'true', '1', 'y', '+'] else 0
    return 0

print("\nConverting binary columns...")
for col in binary_columns:
    if col in df.columns:
        df[col] = df[col].apply(to_binary)
        print(f"{col}: {df[col].sum()} positive values")

# Find significantly upregulated proteins
significant_up = df[
    (df['Log2(CCCP/DMSO)'] > 1) & 
    (df['Adjusted P value (q value; CCCP vs DMSO)'] < 0.05)
]

# Find proteins with ISRIB rescue
isrib_rescue = df[
    (df['Log2(CCCP/DMSO)'] < -1) &  # Decreased by CCCP
    (df['Log2(CCCP+ISRIB/DMSO)'] > df['Log2(CCCP/DMSO)'] + 0.5)  # Improved by ISRIB
]

print(f"\nSignificantly upregulated proteins: {len(significant_up)}")
print(f"Proteins showing ISRIB rescue: {len(isrib_rescue)}")

# Save detailed analysis
with open('proteomics_analysis/results/detailed_analysis.txt', 'w') as f:
    f.write("Detailed Proteomics Analysis\n")
    f.write("==========================\n\n")
    
    # Analysis of upregulated proteins
    f.write("1. Significantly Upregulated Proteins (Log2FC > 1, p < 0.05)\n")
    f.write("-----------------------------------------------------\n")
    f.write(f"Total upregulated proteins: {len(significant_up)}\n\n")
    
    if len(significant_up) > 0:
        for _, row in significant_up.iterrows():
            gene_name = row['Gene name'] if pd.notna(row['Gene name']) else "Unknown"
            location = row['Suborganellar location'] if pd.notna(row['Suborganellar location']) else "Unknown"
            f.write(f"Gene: {gene_name}\n")
            f.write(f"Accession: {row['Accession number']}\n")
            f.write(f"Location: {location}\n")
            f.write(f"Log2FC (CCCP/DMSO): {row['Log2(CCCP/DMSO)']:.3f}\n")
            f.write(f"Log2FC (CCCP+ISRIB/DMSO): {row['Log2(CCCP+ISRIB/DMSO)']:.3f}\n")
            f.write(f"Adjusted p-value: {row['Adjusted P value (q value; CCCP vs DMSO)']:.3e}\n\n")
    
    # Analysis of ISRIB rescue
    f.write("\n2. ISRIB Rescue Analysis\n")
    f.write("----------------------\n")
    f.write(f"Total proteins showing ISRIB rescue: {len(isrib_rescue)}\n")
    mean_rescue = (isrib_rescue['Log2(CCCP+ISRIB/DMSO)'] - isrib_rescue['Log2(CCCP/DMSO)']).mean()
    f.write(f"Mean rescue magnitude: {mean_rescue:.3f}\n")
    
    # Distribution of rescue magnitudes
    rescue_magnitudes = isrib_rescue['Log2(CCCP+ISRIB/DMSO)'] - isrib_rescue['Log2(CCCP/DMSO)']
    f.write(f"\nRescue magnitude statistics:\n")
    f.write(f"25th percentile: {rescue_magnitudes.quantile(0.25):.3f}\n")
    f.write(f"Median: {rescue_magnitudes.median():.3f}\n")
    f.write(f"75th percentile: {rescue_magnitudes.quantile(0.75):.3f}\n")
    f.write(f"Mean: {rescue_magnitudes.mean():.3f}\n")
    f.write(f"Std Dev: {rescue_magnitudes.std():.3f}\n\n")
    
    # Top 10 most rescued proteins
    top_rescue = isrib_rescue.copy()
    top_rescue['Rescue_magnitude'] = top_rescue['Log2(CCCP+ISRIB/DMSO)'] - top_rescue['Log2(CCCP/DMSO)']
    top_rescue = top_rescue.nlargest(10, 'Rescue_magnitude')
    
    f.write("Top 10 Proteins Rescued by ISRIB:\n")
    for _, row in top_rescue.iterrows():
        gene_name = row['Gene name'] if pd.notna(row['Gene name']) else "Unknown"
        location = row['Suborganellar location'] if pd.notna(row['Suborganellar location']) else "Unknown"
        f.write(f"\nGene: {gene_name}\n")
        f.write(f"Accession: {row['Accession number']}\n")
        f.write(f"Location: {location}\n")
        f.write(f"CCCP Log2FC: {row['Log2(CCCP/DMSO)']:.3f}\n")
        f.write(f"CCCP+ISRIB Log2FC: {row['Log2(CCCP+ISRIB/DMSO)']:.3f}\n")
        f.write(f"Rescue magnitude: {row['Rescue_magnitude']:.3f}\n")

    # Category-specific rescue analysis
    f.write("\n3. Category-Specific Rescue Analysis\n")
    f.write("--------------------------------\n")
    for cat in binary_columns:
        if cat in df.columns:
            cat_proteins = df[df[cat] == 1]
            total_cat = len(cat_proteins)
            if total_cat > 0:
                rescued_proteins = cat_proteins[
                    (cat_proteins['Log2(CCCP/DMSO)'] < -1) &
                    (cat_proteins['Log2(CCCP+ISRIB/DMSO)'] > cat_proteins['Log2(CCCP/DMSO)'] + 0.5)
                ]
                
                f.write(f"\n{cat}:\n")
                f.write(f"Total proteins in category: {total_cat}\n")
                f.write(f"Proteins showing rescue: {len(rescued_proteins)}\n")
                
                if len(rescued_proteins) > 0:
                    mean_rescue = (rescued_proteins['Log2(CCCP+ISRIB/DMSO)'] - rescued_proteins['Log2(CCCP/DMSO)']).mean()
                    mean_cccp = rescued_proteins['Log2(CCCP/DMSO)'].mean()
                    mean_cccp_isrib = rescued_proteins['Log2(CCCP+ISRIB/DMSO)'].mean()
                    f.write(f"Mean CCCP Log2FC: {mean_cccp:.3f}\n")
                    f.write(f"Mean CCCP+ISRIB Log2FC: {mean_cccp_isrib:.3f}\n")
                    f.write(f"Mean rescue magnitude: {mean_rescue:.3f}\n")

# Create rescue effect visualization
plt.figure(figsize=(10, 6))
plt.scatter(isrib_rescue['Log2(CCCP/DMSO)'], 
           isrib_rescue['Log2(CCCP+ISRIB/DMSO)'],
           alpha=0.5)
plt.plot([-10, 4], [-10, 4], '--', color='red', alpha=0.5)
plt.xlabel('Log2FC (CCCP/DMSO)')
plt.ylabel('Log2FC (CCCP+ISRIB/DMSO)')
plt.title('ISRIB Rescue Effect')
plt.savefig('proteomics_analysis/plots/isrib_rescue_effect.png', dpi=300, bbox_inches='tight')
plt.close()

# Create rescue magnitude distribution plot
plt.figure(figsize=(10, 6))
rescue_magnitudes = isrib_rescue['Log2(CCCP+ISRIB/DMSO)'] - isrib_rescue['Log2(CCCP/DMSO)']
sns.histplot(data=rescue_magnitudes, kde=True, bins=50)
plt.xlabel('Rescue Magnitude (Log2FC difference)')
plt.ylabel('Count')
plt.title('Distribution of ISRIB Rescue Magnitudes')
plt.savefig('proteomics_analysis/plots/rescue_magnitude_distribution.png', dpi=300, bbox_inches='tight')
plt.close()

print("\nAnalysis complete. Check proteomics_analysis/results/detailed_analysis.txt for results.")
