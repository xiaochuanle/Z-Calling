#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, accuracy_score
import matplotlib.pyplot as plt
import seaborn as sns

motif="TTATT"

### Examplary format of input files (only first 10 rows shown)
"""
Position        Base    ContextKmer     ContextIPD      ContextPW
2031    T       TCTCCCTTTTATCTTTGCCTT   [22, 10, 13, 19, 20, 16, 11, 20, 10, 16, 11, 16, 20, 13, 19, 14, 22, 17, 19, 13, 18]    [15, 20, 22, 19, 15, 10, 10, 12, 24, 19, 9, 12, 11, 13, 13, 24, 25, 14, 12, 22, 37]
4775    A       GTGGTTTTTTATTTGTGGTTT   [20, 11, 18, 19, 21, 21, 16, 15, 14, 11, 18, 15, 28, 19, 22, 21, 15, 18, 13, 15, 17]    [13, 11, 17, 26, 16, 20, 11, 16, 19, 11, 18, 16, 10, 11, 15, 18, 21, 23, 16, 11, 14]
4789    A       GTGGTTTTTTATTTGTGGTTG   [16, 16, 8, 18, 9, 13, 12, 16, 22, 23, 19, 18, 28, 15, 20, 11, 18, 19, 21, 21, 16]      [21, 11, 9, 26, 16, 11, 10, 9, 9, 13, 14, 18, 9, 16, 13, 11, 17, 26, 16, 20, 11]
420     T       TGGTGTTTCTATGTTCTGTTT   [8, 14, 7, 5, 5, 10, 12, 5, 6, 5, 9, 25, 11, 12, 10, 10, 13, 14, 14, 7, 8]      [15, 16, 14, 14, 8, 12, 16, 18, 18, 8, 10, 15, 19, 19, 14, 15, 14, 10, 23, 16, 22]
1655    T       TTTGGGCTTTATTTGTTTTTT   [15, 12, 12, 13, 12, 10, 12, 9, 8, 12, 10, 8, 17, 11, 14, 17, 9, 10, 7, 6, 19]  [15, 18, 20, 23, 17, 16, 18, 11, 10, 15, 16, 20, 9, 11, 18, 18, 19, 19, 11, 10, 13]
2490    A       TTGTCTCTTTATGCCGTTTTC   [14, 15, 12, 14, 20, 14, 9, 7, 10, 14, 9, 11, 15, 14, 11, 13, 16, 18, 17, 7, 17]        [13, 14, 17, 13, 18, 12, 16, 16, 14, 14, 17, 13, 11, 16, 20, 14, 17, 23, 15, 13, 8]
3333    A       CGCTCTCTCTATGTTCTCCTT   [11, 15, 17, 14, 10, 12, 14, 12, 9, 13, 11, 10, 12, 10, 8, 17, 25, 26, 18, 30, 11]      [20, 25, 25, 19, 12, 16, 23, 20, 14, 11, 8, 15, 15, 19, 16, 22, 16, 16, 18, 25, 24]
723     T       CCTTTCGCTTATTTGTTCTGG   [78, 23, 34, 17, 19, 74, 19, 12, 56, 58, 14, 28, 75, 29, 65, 10, 36, 51, 25, 29, 35]    [22, 13, 22, 12, 19, 27, 19, 10, 16, 17, 16, 20, 12, 22, 27, 16, 18, 20, 17, 18, 19]
1807    A       CTCTTCGGGTATCTCTTTTGT   [15, 21, 12, 22, 19, 29, 25, 41, 31, 73, 20, 27, 73, 16, 69, 22, 20, 71, 34, 46, 31]    [24, 18, 10, 8, 14, 18, 15, 20, 14, 27, 12, 19, 13, 11, 18, 14, 18, 15, 20, 31, 22]
"""
file1 = f'm64144_230309_040933.hifi_kinetics.6mA_neg.NT_033779.5.fn3rn3.{motif}.tsv'
file3 = f'm64446e_230310_132726.hifi_kinetics.allZ.NT_033779.5.fn3rn3.{motif}.tsv'

data1 = pd.read_csv(file1, sep='\t')
data3 = pd.read_csv(file3, sep='\t')

# Combine the data into a single DataFrame
data1['Condition'] = 'Canonical A'
data3['Condition'] = 'Z'
combined_data = pd.concat([data1, data3])

# Filter the data for kmers with at least 10 records in both conditions
kmer_counts_data1 = data1['ContextKmer'].value_counts()
kmer_counts_data3 = data3['ContextKmer'].value_counts()
kmers = kmer_counts_data1[kmer_counts_data1 >= 10].index.intersection(
    kmer_counts_data3[kmer_counts_data3 >= 10].index)
combined_data = combined_data[combined_data['ContextKmer'].isin(kmers)].copy()

# Format conversion
combined_data['ContextIPD'] = combined_data['ContextIPD'].apply(eval)
combined_data['ContextPW'] = combined_data['ContextPW'].apply(eval)

def extract_features(row):
    ipd_values = row['ContextIPD']
    pw_values = row['ContextPW']
    features = {}
    for i in range(21):
        features[f'IPD_{i+1}'] = ipd_values[i]
        features[f'PW_{i+1}'] = pw_values[i]
    return pd.Series(features)

# Apply feature extraction
feature_data = combined_data.apply(extract_features, axis=1)
feature_data['Condition'] = combined_data['Condition']
feature_data = feature_data.dropna()

X = feature_data.drop(columns=['Condition'])
y = feature_data['Condition']

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=40)

model = RandomForestClassifier(n_estimators=100, random_state=40)
model.fit(X_train, y_train)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 9

# Predict and evaluate
y_pred = model.predict(X_test)
print(classification_report(y_test, y_pred))
print(f"Accuracy: {accuracy_score(y_test, y_pred)}")

# Analyze feature importance
importances = model.feature_importances_
indices = np.argsort(importances)[::-1]

# Plot feature importances
plt.figure(figsize=(12, 8))
plt.title("Feature Importance (IPD and PW)")
sns.barplot(x=importances[indices], y=[X.columns[i] for i in indices], palette="plasma")
plt.xlabel("Importance")
plt.ylabel("Feature")
plt.show()

