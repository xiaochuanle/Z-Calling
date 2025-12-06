#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

motif = sys.argv[1]

# Load the data
file1 = f'm64144_230309_040933.hifi_kinetics.6mA_neg.NT_033779.5.fn3rn3.{motif}.tsv'
file3 = f'm64446e_230310_132726.hifi_kinetics.allZ.NT_033779.5.fn3rn3.{motif}.tsv'

data1 = pd.read_csv(file1, sep='\t')
data3 = pd.read_csv(file3, sep='\t')

# Count occurrences of 'ContextKmer' in both dataframes
kmer_counts_data1 = data1['ContextKmer'].value_counts()
kmer_counts_data3 = data3['ContextKmer'].value_counts()

# Identify 'ContextKmer' values that occur at least 10 times in both data1 and data3
kmers = kmer_counts_data1[kmer_counts_data1 >= 10].index.intersection(
            kmer_counts_data3[kmer_counts_data3 >= 10].index)

# Filter data1 and data3 to retain only rows with the common 'ContextKmer'
data1 = data1[data1['ContextKmer'].isin(kmers)].copy()
data3 = data3[data3['ContextKmer'].isin(kmers)].copy()

# Function to normalize by the 7th element in the list
#def normalize_signals(series):
#    return series.apply(lambda x: [v / x[6] if x[6] != 0 else np.nan for v in x])

# Normalize ContextIPD and ContextPW values in data1 and data3
data1['ContextIPD'] = data1['ContextIPD'].apply(eval)
data1['ContextPW'] = data1['ContextPW'].apply(eval)
data3['ContextIPD'] = data3['ContextIPD'].apply(eval)
data3['ContextPW'] = data3['ContextPW'].apply(eval)

# Extract IPD and PW values with kmer
def extract_signals(df, condition):
    ipd = df['ContextIPD'].explode().reset_index(drop=True)
    pw = df['ContextPW'].explode().reset_index(drop=True)
    kmer = df['ContextKmer'].repeat(df['ContextIPD'].apply(len)).reset_index(drop=True)
    base_pos = list(range(1, 22)) * (len(ipd) // 21)
    condition = [condition] * len(ipd)
    return pd.DataFrame({'Kmer': kmer, 'BasePosition': base_pos, 'IPD': pd.to_numeric(ipd, errors='coerce'), 'PW': pd.to_numeric(pw, errors='coerce'), 'Condition': condition})

ipd_pw1 = extract_signals(data1, 'Canonical A')
ipd_pw3 = extract_signals(data3, 'Z')

# Combine the data for all conditions
combined_data = pd.concat([ipd_pw1, ipd_pw3])

# Filter out non-numeric values
combined_data = combined_data.dropna()

# Initialize lists to store the results
results = []
p_values_pw = []
p_values_ipd = []

for kmer in kmers:
    kmer_data = combined_data[combined_data['Kmer'] == kmer]
    kmer_result = {'ContextKmer': kmer}
    
    pw_ratios = []
    ipd_ratios = []

    for base_pos in range(1, 22):
        base_data = kmer_data[kmer_data['BasePosition'] == base_pos]
        
        # Split the data into the two conditions
        pw_canonical_a = base_data[base_data['Condition'] == 'Canonical A']['PW']
        pw_z = base_data[base_data['Condition'] == 'Z']['PW']
        ipd_canonical_a = base_data[base_data['Condition'] == 'Canonical A']['IPD']
        ipd_z = base_data[base_data['Condition'] == 'Z']['IPD']

        median_pw_canonical_a = np.median(pw_canonical_a)
        median_pw_z = np.median(pw_z)
        median_ipd_canonical_a = np.median(ipd_canonical_a)
        median_ipd_z = np.median(ipd_z)
        
        pw_ratios.append(median_pw_z / median_pw_canonical_a if median_pw_canonical_a != 0 else np.nan)
        ipd_ratios.append(median_ipd_z / median_ipd_canonical_a if median_ipd_canonical_a != 0 else np.nan)

    pw_norm_factor = np.median(pw_ratios)
    ipd_norm_factor = np.median(ipd_ratios)

    for base_pos in range(1, 22):
        base_data = kmer_data[kmer_data['BasePosition'] == base_pos]
        
        # Split the data into the two conditions
        pw_canonical_a = base_data[base_data['Condition'] == 'Canonical A']['PW']
        pw_z = base_data[base_data['Condition'] == 'Z']['PW'] / pw_norm_factor
        ipd_canonical_a = base_data[base_data['Condition'] == 'Canonical A']['IPD']
        ipd_z = base_data[base_data['Condition'] == 'Z']['IPD'] / ipd_norm_factor

        median_pw_canonical_a = np.median(pw_canonical_a)
        median_pw_z = np.median(pw_z)
        median_ipd_canonical_a = np.median(ipd_canonical_a)
        median_ipd_z = np.median(ipd_z)
        
        pw_ratio = median_pw_z / median_pw_canonical_a if median_pw_canonical_a != 0 else np.nan
        ipd_ratio = median_ipd_z / median_ipd_canonical_a if median_ipd_canonical_a != 0 else np.nan

        # Perform two-tailed t-tests for PW and IPD
        t_pw, p_pw = ttest_ind(pw_canonical_a, pw_z, nan_policy='omit')
        t_ipd, p_ipd = ttest_ind(ipd_canonical_a, ipd_z, nan_policy='omit')
        
        p_values_pw.append(p_pw)
        p_values_ipd.append(p_ipd)
        
        # Calculate median ratios for PW and IPD
        kmer_result[f'PW_pvalue_BasePosition_{base_pos}'] = p_pw
        kmer_result[f'IPD_pvalue_BasePosition_{base_pos}'] = p_ipd
        kmer_result[f'PW_Ratio_BasePosition_{base_pos}'] = pw_ratio
        kmer_result[f'IPD_Ratio_BasePosition_{base_pos}'] = ipd_ratio
    
    results.append(kmer_result)

# BH correction for p-values
#p_values_pw_corrected = multipletests(p_values_pw, method='fdr_bh')[1]
#p_values_ipd_corrected = multipletests(p_values_ipd, method='fdr_bh')[1]

import math
p_values_pw_corrected = multipletests([1.0 if math.isnan(x) else x for x in p_values_pw], method='fdr_bh')[1]
p_values_ipd_corrected = multipletests([1.0 if math.isnan(x) else x for x in p_values_ipd], method='fdr_bh')[1]

# Insert corrected p-values back into the results
for i, kmer_result in enumerate(results):
    for base_pos in range(1, 22):
        kmer_result[f'PW_pvalue_BasePosition_{base_pos}'] = p_values_pw_corrected[i * 21 + (base_pos - 1)]
        kmer_result[f'IPD_pvalue_BasePosition_{base_pos}'] = p_values_ipd_corrected[i * 21 + (base_pos - 1)]

# Convert results to DataFrame and save
final_df = pd.DataFrame(results)
output_file = f"{motif}_kmer_analysis_results.tsv"
final_df.to_csv(output_file, sep='\t', index=False)
