#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
from random import sample

# Load the data
file_A = 'm64144_230309_040933.hifi_kinetics.6mA_neg.NT_033779.5.fn3rn3.tsv'  # Replace with your file path for Condition A
file_Z = 'm64446e_230310_132726.hifi_kinetics.allZ.NT_033779.5.fn3rn3.tsv'  # Replace with your file path for Condition Z

df_A = pd.read_csv(file_A, sep='\t')
df_Z = pd.read_csv(file_Z, sep='\t')

# Add the 'motif' column (the center 3 bases from the ContextKmer)
df_A['motif'] = df_A['ContextKmer'].apply(lambda x: x[9:12])
df_Z['motif'] = df_Z['ContextKmer'].apply(lambda x: x[9:12])

# Randomly downsample 10000 rows per unique motif value in each dataframe
def downsample(df, n=10000):
    return df.groupby('motif').apply(lambda x: x.sample(n=min(len(x), n), random_state=42)).reset_index(drop=True)

df_A_downsampled = downsample(df_A)
df_Z_downsampled = downsample(df_Z)

# Base encoding: A=1, C=2, G=3, T=4
base_encoding = {'A': 1, 'C': 2, 'G': 3, 'T': 4}

# Encode the ContextKmer to a sequence of integers
df_A_downsampled['EncodedKmer'] = df_A_downsampled['ContextKmer'].apply(lambda x: [base_encoding[base] for base in x])
df_Z_downsampled['EncodedKmer'] = df_Z_downsampled['ContextKmer'].apply(lambda x: [base_encoding[base] for base in x])

df_A_downsampled['Condition'] = "A"
df_Z_downsampled['Condition'] = "Z"

# Concatenate the downsampled dataframes
df_combined = pd.concat([df_A_downsampled, df_Z_downsampled])

# Prepare data for PCA and t-SNE
X = np.concatenate([df_combined['ContextIPD'].apply(lambda x: eval(x)[5:16]).values.tolist(),
                    df_combined['ContextPW'].apply(lambda x: eval(x)[5:16]).values.tolist(),
                    df_combined['EncodedKmer'].apply(lambda x: x[5:16]).values.tolist()], axis=1)

motif = ["TAT", "CAC", "GAG", "TAC", "TAG", "CAT", "CAG", "GAT", "GAC"]
df_combined_motif = df_combined.loc[df_combined['motif'].isin(motif)]

# Perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(X)
df_combined_motif['PCA1'] = pca_result[:, 0]
df_combined_motif['PCA2'] = pca_result[:, 1]

# Perform t-SNE
tsne = TSNE(n_components=2, random_state=42)
tsne_result = tsne.fit_transform(X)
df_combined_motif['TSNE1'] = tsne_result[:, 0]
df_combined_motif['TSNE2'] = tsne_result[:, 1]

from matplotlib.backends.backend_pdf import PdfPages

# Set font parameters (if not already set in your environment)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 7

# Create a PdfPages object to save the figures
with PdfPages('tsne_pca_plots.pdf') as pdf:

    # First plot: PCA1 vs PCA2, colored by Condition
    g = sns.jointplot(x='PCA1', y='PCA2', hue='Condition', data=downsample(df_combined_motif, 1000), palette='Set1', joint_kws={'s':10}, marginal_kws={'alpha':0.05, 'linewidth':0.3}, alpha=0.6)
    g.fig.set_size_inches(4, 3)  # Set figure size to 8x6 inches
    g.ax_joint.legend(loc='upper left', bbox_to_anchor=(0.85, 1), frameon=False, framealpha=0.3, handletextpad=0.05, labelspacing=0.4, borderaxespad=0.3)
    g.fig.suptitle('PCA: Colored by Condition')
    g.fig.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to prevent clipping of the title and axis labels
    pdf.savefig(g.fig, bbox_inches='tight')  # Save figure with tight bounding box
    plt.close(g.fig)  # Close the figure after saving

    # Second plot: PCA1 vs PCA2, colored by Motif
    g = sns.jointplot(x='PCA1', y='PCA2', hue='motif', data=downsample(df_combined_motif, 1000), palette='Set2', joint_kws={'s':10}, marginal_kws={'alpha':0.05, 'linewidth':0.3}, alpha=0.6)
    g.fig.set_size_inches(4, 3)  # Set figure size to 8x6 inches
    g.ax_joint.legend(loc='upper left', bbox_to_anchor=(0.8, 1), frameon=False, framealpha=0.3, handletextpad=0.05, labelspacing=0.4, borderaxespad=0.3)
    g.fig.suptitle('PCA: Colored by Motif')
    g.fig.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to prevent clipping of the title and axis labels
    pdf.savefig(g.fig, bbox_inches='tight')  # Save figure with tight bounding box
    plt.close(g.fig)  # Close the figure after saving

    # Third plot: TSNE1 vs TSNE2, colored by Condition
    g = sns.jointplot(x='TSNE1', y='TSNE2', hue='Condition', data=downsample(df_combined_motif, 1000), palette='Set1', joint_kws={'s':10}, marginal_kws={'alpha':0.05, 'linewidth':0.3}, alpha=0.6)
    g.fig.set_size_inches(4, 3)  # Set figure size to 8x6 inches
    g.ax_joint.legend(loc='upper left', bbox_to_anchor=(0.85, 1), frameon=False, framealpha=0.3, handletextpad=0.05, labelspacing=0.4, borderaxespad=0.3)
    g.fig.suptitle('t-SNE: Colored by Condition')
    g.fig.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to prevent clipping of the title and axis labels
    pdf.savefig(g.fig, bbox_inches='tight')  # Save figure with tight bounding box
    plt.close(g.fig)  # Close the figure after saving

    # Fourth plot: TSNE1 vs TSNE2, colored by Motif
    g = sns.jointplot(x='TSNE1', y='TSNE2', hue='motif', data=downsample(df_combined_motif, 1000), palette='Set2', joint_kws={'s':10}, marginal_kws={'alpha':0.05, 'linewidth':0.3}, alpha=0.6)
    g.fig.set_size_inches(4, 3)  # Set figure size to 8x6 inches
    g.ax_joint.legend(loc='upper left', bbox_to_anchor=(0.8, 1), frameon=False, framealpha=0.3, handletextpad=0.05, labelspacing=0.4, borderaxespad=0.3)
    g.fig.suptitle('t-SNE: Colored by Motif')
    g.fig.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to prevent clipping of the title and axis labels
    pdf.savefig(g.fig, bbox_inches='tight')  # Save figure with tight bounding box
    plt.close(g.fig)  # Close the figure after saving