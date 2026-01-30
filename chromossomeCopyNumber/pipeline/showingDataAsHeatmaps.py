import pandas as pd
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt

# List all bedgraph files in the 'output' directory
file_list = glob.glob('output/*mapbed_output.bedgraph')

# Check if any files were found
if not file_list:
    print("No BEDGRAPH files found in the 'output' directory.")
    exit(1)

# Prepare a DataFrame to hold all results
results_df = pd.DataFrame()

for file_path in file_list:
    # Read the bedgraph file with 4 columns: chrom, start, end, score
    df = pd.read_csv(file_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])

    # Skip files with empty or NaN 'score' column
    if df['score'].isna().all():
        print(f"Warning: File {file_path} contains only NaN values in the 'score' column. Skipping this file.")
        continue

    # Drop rows where 'score' is NaN to avoid issues with calculations
    df = df.dropna(subset=['score'])

    # Calculate CCN
    median_score = df['score'].median()
    if median_score == 0:
        print(f"Warning: Median 'score' is 0 in the file {file_path}. Skipping this file.")
        continue

    df['CCN'] = (df['score'] / median_score) * 2

    # Aggregate the CCN value by chromosome (scaffold)
    agg_df = df.groupby('chrom')['CCN'].median().reset_index()

    # Extract the filename (without extension) to use as the column header
    column_header = os.path.splitext(os.path.basename(file_path))[0]

    # If it's the first file, initialize the results DataFrame with chromosome names
    if results_df.empty:
        results_df = agg_df.rename(columns={'CCN': column_header})
    else:
        # For subsequent files, add the results as a new column
        results_df = results_df.merge(agg_df, on='chrom', how='outer').rename(columns={'CCN': column_header})

# Save the aggregated results to a TSV file
os.makedirs('output', exist_ok=True)
results_df.to_csv('output/bedgraph_aggregated_results.tsv', sep='\t', index=False)

# Create a heatmap of the results
plt.figure(figsize=(14, 10))
heatmap = sns.heatmap(results_df.set_index('chrom').transpose(), cmap='viridis', vmin=0, vmax=4, cbar_kws={'label': 'CCN'})
heatmap.set_title('Chromosome CCN Heatmap')
heatmap.set_xlabel('Chromosome')
heatmap.set_ylabel('Sample')
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()  # Adjust spacing

# Save the heatmap as a PNG image
plt.savefig('output/CCN_heatmap.png')

print("Processing complete. Results saved to 'output/bedgraph_aggregated_results.tsv' and 'output/CCN_heatmap.png'.")
