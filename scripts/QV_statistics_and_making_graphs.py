#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# Assuming you have two DataFrames df_canonical_a and df_z for conditions "Canonical A" and "Z"
# Each DataFrame has 'np' and 'QV' columns

# Add a column to each DataFrame to indicate the condition
revio_z = pd.read_table('/data1/FRUITFLY/REVIO_S2ZTP/REVIO.S2ZTP.10000.tsv')
revio_a = pd.read_table('/data1/FRUITFLY/REVIO_S2NEG/REVIO.S2NEG.10000.tsv')

revio_a = revio_a.loc[(revio_a.QV < 100)].copy()

revio_z = revio_z.loc[(revio_z.QV < 100)].copy()

# Add a column to each DataFrame to indicate the condition
revio_a['Condition'] = 'Canonical A'
revio_z['Condition'] = 'Z'
# Combine the two DataFrames
df_combined = pd.concat([revio_a, revio_z])

df_combined = df_combined.reset_index(drop=True)

df_combined['Pacbio'] = 'Revio'

# Add a column to each DataFrame to indicate the condition
seq2_z = pd.read_table('/data1/FRUITFLY/SEQ2_S2ZTP/SEQ2.S2ZTP.10000.tsv')
seq2_a = pd.read_table('/data1/FRUITFLY/SEQ2_S2NEG/SEQ2.S2NEG.10000.tsv')

seq2_a = seq2_a.loc[(seq2_a.QV < 100)].copy()
seq2_z = seq2_z.loc[(seq2_z.QV < 100)].copy()

# Add a column to each DataFrame to indicate the condition
seq2_a['Condition'] = 'Canonical A'
seq2_z['Condition'] = 'Z'
# Combine the two DataFrames
df2_combined = pd.concat([seq2_a, seq2_z])
df2_combined = df2_combined.reset_index(drop=True)

df2_combined['Pacbio'] = 'Sequel_2'

df = pd.concat([df_combined,df2_combined]).reset_index(drop=True)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 12

pdf_output_file = f"/data1/FRUITFLY/Z_A_QV_lineplots.pdf"

plt.figure(figsize=(6, 2))
sns.relplot(data=df, x='NP', y='QV', col='Pacbio', hue='Condition', style='Condition', errorbar='sd', markers=True, dashes = True, kind='line')

# Set the x-axis to use integer ticks only
plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))

# Show plot
plt.savefig(pdf_output_file)
plt.close()

# Separate the data by Pacbio platform
platforms = df['Pacbio'].unique()
results = []

for platform in platforms:
    platform_data = df[df['Pacbio'] == platform]
    p_values = []
    np_values = []
    mean_qvs_a = []
    mean_qvs_z = []
    std_qvs_a = []
    std_qvs_z = []

    # Perform Mann-Whitney U tests and calculate statistics for each NP value
    for np_value in platform_data['NP'].unique():
        subset = platform_data[platform_data['NP'] == np_value]
        group_a = subset[subset['Condition'] == 'Canonical A']['QV']
        group_z = subset[subset['Condition'] == 'Z']['QV']

        # Perform the Mann-Whitney U test
        u_stat, p_value = mannwhitneyu(group_a, group_z, alternative='two-sided')
        p_values.append(p_value)
        np_values.append(np_value)

        # Calculate mean and std deviation for QV
        mean_qvs_a.append(group_a.mean())
        mean_qvs_z.append(group_z.mean())
        std_qvs_a.append(group_a.std())
        std_qvs_z.append(group_z.std())

    # Apply BH correction to the p-values
    adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]

    # Store the results
    platform_results = pd.DataFrame({
        'NP': np_values,
        f'{platform}_adjusted_p_value': adjusted_p_values,
        f'{platform}_mean_QV_Canonical_A': mean_qvs_a,
        f'{platform}_std_QV_Canonical_A': std_qvs_a,
        f'{platform}_mean_QV_Z': mean_qvs_z,
        f'{platform}_std_QV_Z': std_qvs_z
    })

    results.append(platform_results)

# Merge results for both platforms
final_results = pd.concat(results, axis=1)

# Remove duplicate NP columns and align the results correctly
final_results = final_results.loc[:, ~final_results.columns.duplicated()]

# Output the final DataFrame with adjusted p-values and statistics
print(final_results)

# Save the final results to a TSV file
final_results.to_csv('/data1/FRUITFLY/adjusted_p_values_and_stats_by_platform_mannwhitneyu.tsv', sep='\t', index=False)