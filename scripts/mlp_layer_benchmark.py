#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score, recall_score, confusion_matrix
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import sys
import argparse
import os
import warnings

# --- CONFIGURATION & STYLE ---
# Suppress Warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

# Set plotting style for Communications Biology (Clean, minimal)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
sns.set_context("paper", font_scale=1.4)
sns.set_style("ticks")

# --- ARGUMENT PARSING ---
parser = argparse.ArgumentParser(description="Benchmark MLP Hidden Layers for Z-base detection.")
parser.add_argument("file_a", help="Input file for Canonical A (TSV)")
parser.add_argument("file_z", help="Input file for Z-base (TSV)")
parser.add_argument("motif", help="Motif to filter by (e.g., TTATT)")
parser.add_argument("-s", "--sample_size", type=int, default=100000, 
                    help="Number of samples per class for training. Default: 100000")
parser.add_argument("-k", "--kmer_length", type=int, default=21, 
                    help="Length of center k-mer to use. Default: 21")
parser.add_argument("-b", "--batch_size", type=int, default=2048,
                    help="Batch size for training and evaluation. Default: 2048")

args = parser.parse_args()

FILE_A = args.file_a
FILE_Z = args.file_z
MOTIF = args.motif
TRAIN_SAMPLE_SIZE = args.sample_size
KMER_LENGTH = args.kmer_length
BATCH_SIZE = args.batch_size
INPUT_SEQ_LEN = 21 
SLICE_START = (INPUT_SEQ_LEN - KMER_LENGTH) // 2
SLICE_END = SLICE_START + KMER_LENGTH

# --- DATA PROCESSING HELPERS ---
def parse_kinetic_string_and_slice(s):
    try:
        if isinstance(s, str):
            full_array = np.array(eval(s))
        else:
            full_array = np.array(s)
        if len(full_array) != INPUT_SEQ_LEN:
            if len(full_array) == 0: return np.zeros(KMER_LENGTH)
        return full_array[SLICE_START:SLICE_END]
    except:
        return np.zeros(KMER_LENGTH)

def one_hot_encode_sequence_slice(seq):
    sliced_seq = seq[SLICE_START:SLICE_END]
    mapping = {'A': [1,0,0,0], 'C': [0,1,0,0], 'G': [0,0,1,0], 'T': [0,0,0,1], 'N': [0,0,0,0]}
    encoded = []
    for base in sliced_seq:
        encoded.extend(mapping.get(base, [0,0,0,0]))
    return np.array(encoded)

def process_single_file(filepath, label, motif_filter=None):
    print(f"Loading {filepath}...")
    try:
        df = pd.read_csv(filepath, sep='\t')
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None, None
    
    if motif_filter and 'ContextKmer' in df.columns:
        df = df[df['ContextKmer'].str.contains(motif_filter)]
    
    if df.empty: return None, None

    ipd_data = np.stack(df['ContextIPD'].apply(parse_kinetic_string_and_slice).values)
    pw_data = np.stack(df['ContextPW'].apply(parse_kinetic_string_and_slice).values)
    seq_data = np.stack(df['ContextKmer'].apply(one_hot_encode_sequence_slice).values)
    
    X_part = np.hstack([seq_data, ipd_data, pw_data])
    y_part = np.full(X_part.shape[0], label)
    return X_part, y_part

def load_data():
    X_a, y_a = process_single_file(FILE_A, 0, MOTIF)
    X_z, y_z = process_single_file(FILE_Z, 1, MOTIF)
    if X_a is None or X_z is None: sys.exit("Error: Empty data.")
    return np.vstack([X_a, X_z]), np.concatenate([y_a, y_z])

def sample_balanced(X, y, size, seed):
    rng = np.random.RandomState(seed)
    idx_0 = np.where(y == 0)[0]
    idx_1 = np.where(y == 1)[0]
    n_0 = min(len(idx_0), size)
    n_1 = min(len(idx_1), size)
    sub_0 = rng.choice(idx_0, n_0, replace=False)
    sub_1 = rng.choice(idx_1, n_1, replace=False)
    indices = np.concatenate([sub_0, sub_1])
    rng.shuffle(indices)
    return X[indices], y[indices]

# --- DYNAMIC MODEL ARCHITECTURE ---
class DynamicMLP(nn.Module):
    def __init__(self, input_dim, num_layers, hidden_dim=128):
        super(DynamicMLP, self).__init__()
        layers = []
        in_d = input_dim
        
        for _ in range(num_layers):
            layers.append(nn.Linear(in_d, hidden_dim))
            layers.append(nn.BatchNorm1d(hidden_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(0.2))
            in_d = hidden_dim
            
        layers.append(nn.Linear(in_d, 1)) # Output Logits
        self.network = nn.Sequential(*layers)

    def forward(self, x):
        return self.network(x)

# --- TRAINING ENGINE ---
def train_and_evaluate(X_train, y_train, X_test, y_test, num_layers, device):
    input_dim = X_train.shape[1]
    model = DynamicMLP(input_dim, num_layers).to(device)
    
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    # Train DataLoader
    train_ds = TensorDataset(torch.FloatTensor(X_train), torch.FloatTensor(y_train).unsqueeze(1))
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    
    # Test DataLoader (For batched inference to save VRAM)
    test_ds = TensorDataset(torch.FloatTensor(X_test))
    test_loader = DataLoader(test_ds, batch_size=BATCH_SIZE, shuffle=False)
    
    # Simple Training Loop with Early Stopping
    model.train()
    max_epochs = 50 
    best_loss = float('inf')
    
    for epoch in range(max_epochs):
        epoch_loss = 0
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            outputs = model(xb)
            loss = criterion(outputs, yb)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
        
        if epoch_loss < best_loss:
            best_loss = epoch_loss
    
    # Evaluation (Batched Inference)
    model.eval()
    all_logits = []
    
    with torch.no_grad():
        for (xb,) in test_loader:
            xb = xb.to(device)
            logits = model(xb)
            all_logits.append(logits.cpu().numpy())
            
    # Concatenate all batches
    logits_flat = np.vstack(all_logits).flatten()
    probs = 1 / (1 + np.exp(-logits_flat)) # Sigmoid
    preds = (probs > 0.5).astype(int)
    
    auc = roc_auc_score(y_test, probs)
    acc = accuracy_score(y_test, preds)
    tn, fp, fn, tp = confusion_matrix(y_test, preds).ravel()
    sensitivity = tp / (tp + fn) # Recall
    specificity = tn / (tn + fp)
    
    return auc, acc, sensitivity, specificity

# --- MAIN BENCHMARK LOOP ---
def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    print("Loading and merging data...")
    X_full, y_full = load_data()
    
    scaler = StandardScaler()
    X_scaled_full = scaler.fit_transform(X_full)
    
    layer_options = [1, 2, 3, 4, 5]
    replicate_seeds = [42, 101, 999]
    results = []
    
    print(f"\nStarting Benchmark: Layers {layer_options}, Replicates {len(replicate_seeds)}")
    print("-" * 60)
    print(f"{'Layers':<10} | {'Rep':<5} | {'AUC':<10} | {'Acc':<10}")
    print("-" * 60)
    
    for n_layers in layer_options:
        for seed in replicate_seeds:
            X_train, X_test, y_train, y_test = train_test_split(
                X_scaled_full, y_full, test_size=0.3, random_state=seed, stratify=y_full
            )
            
            X_sub, y_sub = sample_balanced(X_train, y_train, TRAIN_SAMPLE_SIZE, seed)
            
            auc, acc, sens, spec = train_and_evaluate(
                X_sub, y_sub, X_test, y_test, n_layers, device
            )
            
            print(f"{n_layers:<10} | {seed:<5} | {auc:.4f}     | {acc:.4f}")
            
            results.append({
                'Hidden_Layers': n_layers,
                'Replicate_Seed': seed,
                'AUC': auc,
                'Accuracy': acc,
                'Sensitivity': sens,
                'Specificity': spec
            })
            
    # --- SAVE RESULTS ---
    df_results = pd.DataFrame(results)
    out_prefix = f"{MOTIF}_k{KMER_LENGTH}_s{TRAIN_SAMPLE_SIZE}_mlp_benchmark"
    df_results.to_csv(f"{out_prefix}_metrics.tsv", sep='\t', index=False)
    
    # --- PLOTTING FOR MANUSCRIPT ---
    # Reshape (Melt) Data for Hue Plotting
    df_melted = df_results.melt(
        id_vars=['Hidden_Layers', 'Replicate_Seed'], 
        value_vars=['AUC', 'Accuracy'],
        var_name='Metric', 
        value_name='Score'
    )
    
    # Palette: Blue (AUC), Vermilion (Accuracy)
    palette = {'AUC': '#0072B2', 'Accuracy': '#D55E00'} 

    plt.figure(figsize=(8, 6))
    
    # 1. Pointplot: Draws Mean + SD Error Bars + Connecting Lines
    # Use 'dodge' to slightly separate AUC and Accuracy visually
    sns.pointplot(
        data=df_melted, x='Hidden_Layers', y='Score', hue='Metric',
        palette=palette, errorbar='sd', capsize=0.1, dodge=0.2,
        markers=['o', 's'], linestyles=['-', '--'], scale=0.8
    )

    # 2. Stripplot: Overlays the raw data points
    # Must use same 'dodge' value to align points with the error bars
    sns.stripplot(
        data=df_melted, x='Hidden_Layers', y='Score', hue='Metric',
        palette=palette, dodge=0.2, size=5, alpha=0.6, jitter=True,
        legend=False # Disable legend for stripplot to avoid duplicates
    )
    
    plt.title(f'Effect of MLP Depth on Z-Calling Performance\n({MOTIF}, {KMER_LENGTH}-mer)', fontsize=14, fontweight='bold')
    plt.xlabel('Number of Hidden Layers', fontsize=12, fontweight='bold')
    plt.ylabel('Score (Mean ± SD)', fontsize=12, fontweight='bold')
    
    # Adjust Y-limits
    min_val = df_melted['Score'].min()
    plt.ylim(max(0.75, min_val - 0.02), 1.0)
    
    sns.despine()
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.legend(frameon=False, fontsize=12, loc='lower right', title=None)
    plt.tight_layout()
    
    plt.savefig(f"{out_prefix}_performance.pdf", format='pdf', dpi=300)
    print(f"\nBenchmark Complete. Files saved:\n1. {out_prefix}_metrics.tsv\n2. {out_prefix}_performance.pdf")

if __name__ == "__main__":
    main()
