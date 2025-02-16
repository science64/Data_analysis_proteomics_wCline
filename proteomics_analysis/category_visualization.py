import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set the default style
plt.style.use('default')

# Read the data
df = pd.read_excel("1-s2.0-S1097276521009540-mmc6.xlsx", sheet_name='Uptake')

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

# Define categories and their display names
categories = {
    'TOM complex': 'TOM Complex',
    'Respiratory complex I': 'Complex I',
    'Mitochondrial small ribosomal subunit': 'Mito Ribosomes',
    'Mitochondrial genome-encoded proteins': 'Mt-encoded'
}

print("Processing categories...")

# Convert binary columns
for cat in categories.keys():
    if cat in df.columns:
        df[cat] = df[cat].apply(to_binary)
        print(f"{cat}: {df[cat].sum()} positive values")

# Prepare data for visualization
category_stats = []
for cat, display_name in categories.items():
    if cat in df.columns:
        cat_proteins = df[df[cat] == 1]
        total = len(cat_proteins)
        if total > 0:
            # Calculate statistics
            rescued = len(cat_proteins[
                (cat_proteins['Log2(CCCP/DMSO)'] < -1) &
                (cat_proteins['Log2(CCCP+ISRIB/DMSO)'] > cat_proteins['Log2(CCCP/DMSO)'] + 0.5)
            ])
            mean_cccp = cat_proteins['Log2(CCCP/DMSO)'].mean()
            mean_cccp_isrib = cat_proteins['Log2(CCCP+ISRIB/DMSO)'].mean()
            rescue_magnitude = mean_cccp_isrib - mean_cccp
            
            print(f"\n{display_name}:")
            print(f"Total proteins: {total}")
            print(f"Rescued: {rescued}")
            print(f"Mean CCCP Log2FC: {mean_cccp:.3f}")
            print(f"Mean CCCP+ISRIB Log2FC: {mean_cccp_isrib:.3f}")
            
            category_stats.append({
                'Category': display_name,
                'Total': total,
                'Rescued': rescued,
                'Rescue_Percentage': (rescued/total)*100,
                'Mean_CCCP': mean_cccp,
                'Mean_CCCP_ISRIB': mean_cccp_isrib,
                'Rescue_Magnitude': rescue_magnitude
            })

# Convert to DataFrame
stats_df = pd.DataFrame(category_stats)
print("\nGenerated statistics DataFrame:")
print(stats_df)

print("\nCreating visualizations...")

# Create visualization 1: Rescue Percentages
plt.figure(figsize=(10, 6))
colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']
bars = plt.bar(range(len(stats_df)), stats_df['Rescue_Percentage'], color=colors)
plt.axhline(y=50, color='r', linestyle='--', alpha=0.3, label='50% threshold')
plt.title('Percentage of Proteins Rescued by ISRIB in Each Category', pad=20)
plt.ylabel('Percentage Rescued (%)')
plt.xticks(range(len(stats_df)), stats_df['Category'], rotation=45, ha='right')

# Add value labels on bars
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height,
             f'{height:.1f}%',
             ha='center', va='bottom')

plt.legend()
plt.tight_layout()
plt.savefig('proteomics_analysis/plots/category_rescue_percentages.png', dpi=300, bbox_inches='tight')
plt.close()

# Create visualization 2: Treatment Effects
plt.figure(figsize=(10, 6))
x = np.arange(len(stats_df))
width = 0.35

plt.bar(x - width/2, stats_df['Mean_CCCP'], width, label='CCCP', color='#ff9999', alpha=0.8)
plt.bar(x + width/2, stats_df['Mean_CCCP_ISRIB'], width, label='CCCP+ISRIB', color='#66b3ff', alpha=0.8)

plt.xlabel('Protein Categories')
plt.ylabel('Mean Log2 Fold Change')
plt.title('CCCP vs CCCP+ISRIB Effects by Category', pad=20)
plt.xticks(x, stats_df['Category'], rotation=45, ha='right')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('proteomics_analysis/plots/category_treatment_effects.png', dpi=300, bbox_inches='tight')
plt.close()

# Create visualization 3: Rescue Magnitudes
plt.figure(figsize=(10, 6))
colors = ['#ff9999' if x < 0 else '#66b3ff' for x in stats_df['Rescue_Magnitude']]
bars = plt.bar(range(len(stats_df)), stats_df['Rescue_Magnitude'], color=colors)
plt.title('ISRIB Rescue Magnitude by Category', pad=20)
plt.ylabel('Mean Rescue Magnitude (Log2FC)')
plt.xticks(range(len(stats_df)), stats_df['Category'], rotation=45, ha='right')
plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

# Add value labels on bars
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height,
             f'{height:.2f}',
             ha='center', va='bottom' if height > 0 else 'top')

plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('proteomics_analysis/plots/category_rescue_magnitudes.png', dpi=300, bbox_inches='tight')
plt.close()

# Create visualization 4: Combined metrics heatmap
metrics_df = stats_df.copy()
metrics_df['Relative_Impact'] = -metrics_df['Mean_CCCP']  # Negative because more negative means bigger impact
metrics_df['Recovery_Level'] = metrics_df['Mean_CCCP_ISRIB']
metrics_df['Rescue_Efficiency'] = metrics_df['Rescue_Magnitude']

heatmap_data = metrics_df[['Category', 'Relative_Impact', 'Recovery_Level', 'Rescue_Efficiency']].set_index('Category')
heatmap_data = (heatmap_data - heatmap_data.mean()) / heatmap_data.std()  # Z-score normalization

plt.figure(figsize=(10, 6))
sns.heatmap(heatmap_data, annot=True, cmap='RdBu_r', center=0, fmt='.2f',
            cbar_kws={'label': 'Z-score'})
plt.title('Normalized Metrics by Category', pad=20)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('proteomics_analysis/plots/category_metrics_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()

# Save numerical results
with open('proteomics_analysis/results/category_analysis_summary.txt', 'w') as f:
    f.write("Category Analysis Summary\n")
    f.write("=======================\n\n")
    
    for _, row in stats_df.iterrows():
        f.write(f"{row['Category']}:\n")
        f.write(f"- Total proteins: {row['Total']}\n")
        f.write(f"- Proteins rescued: {row['Rescued']} ({row['Rescue_Percentage']:.1f}%)\n")
        f.write(f"- Mean CCCP effect: {row['Mean_CCCP']:.3f}\n")
        f.write(f"- Mean CCCP+ISRIB effect: {row['Mean_CCCP_ISRIB']:.3f}\n")
        f.write(f"- Rescue magnitude: {row['Rescue_Magnitude']:.3f}\n\n")

print("\nCategory visualizations complete. Check proteomics_analysis/plots/ directory for the outputs.")
