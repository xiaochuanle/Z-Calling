#!/usr/bin/env python3

import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import sys
import argparse
import os
import warnings

# Suppress Warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

# --- CONFIGURATION ---
parser = argparse.ArgumentParser(description="Run SHAP analysis on PyTorch MLP model for Z-base detection.")
parser.add_argument("file_a", help="Input file for Canonical A (TSV)")
parser.add_argument("file_z", help="Input file for Z-base (TSV)")
parser.add_argument("motif", help="Motif to filter by (e.g., TTATT)")
parser.add_argument("-s", "--sample_size", type=int, default=10000, 
                    help="Number of samples per class for training. Default: 10000")
parser.add_argument("-k", "--kmer_length", type=int, default=21, 
                    help="Length of center k-mer to use (odd number, <= 21). Default: 21")

args = parser.parse_args()

FILE_A = args.file_a
FILE_Z = args.file_z
MOTIF = args.motif
TRAIN_SAMPLE_SIZE = args.sample_size
KMER_LENGTH = args.kmer_length
RANDOM_SEED = 42

INPUT_SEQ_LEN = 21 

if KMER_LENGTH % 2 == 0:
    print(f"Error: K-mer length must be odd. Got {KMER_LENGTH}.")
    sys.exit(1)
if KMER_LENGTH > INPUT_SEQ_LEN:
    print(f"Error: K-mer length ({KMER_LENGTH}) cannot be larger than input length ({INPUT_SEQ_LEN}).")
    sys.exit(1)

SLICE_START = (INPUT_SEQ_LEN - KMER_LENGTH) // 2
SLICE_END = SLICE_START + KMER_LENGTH

print(f"Configuration: Motif={MOTIF}, K-mer Length={KMER_LENGTH}, Train Size={TRAIN_SAMPLE_SIZE}")
print(f"  - Slicing input arrays from index {SLICE_START} to {SLICE_END}")

# --- DEVICE CONFIGURATION (CUDA SUPPORT) ---
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"  - PyTorch computing device: {device}")

# --- DATA PROCESSING ---

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
    mapping = {
        'A': [1, 0, 0, 0],
        'C': [0, 1, 0, 0],
        'G': [0, 0, 1, 0],
        'T': [0, 0, 0, 1],
        'N': [0, 0, 0, 0]
    }
    encoded = []
    for base in sliced_seq:
        encoded.extend(mapping.get(base, [0, 0, 0, 0]))
    return np.array(encoded)

def process_single_file(filepath, label, motif_filter=None):
    print(f"Loading {filepath}...")
    try:
        df = pd.read_csv(filepath, sep='\t')
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None, None, None, None

    if motif_filter and 'ContextKmer' in df.columns:
        original_len = len(df)
        df = df[df['ContextKmer'].str.contains(motif_filter)]
        print(f"  - Filtered {original_len} -> {len(df)} rows with motif '{motif_filter}'")

    if df.empty:
        return None, None, None, None

    ipd_data = np.stack(df['ContextIPD'].apply(parse_kinetic_string_and_slice).values)
    pw_data = np.stack(df['ContextPW'].apply(parse_kinetic_string_and_slice).values)
    
    if 'ContextKmer' not in df.columns:
        raise ValueError(f"File {filepath} missing 'ContextKmer' column required for sequence features.")
        
    seq_data = np.stack(df['ContextKmer'].apply(one_hot_encode_sequence_slice).values)

    X_part = np.hstack([seq_data, ipd_data, pw_data])
    y_part = np.full(X_part.shape[0], label)
    
    return X_part, y_part, KMER_LENGTH

def load_and_merge_files(file_a, file_z):
    X_a, y_a, len_a = process_single_file(file_a, label=0, motif_filter=MOTIF)
    X_z, y_z, len_z = process_single_file(file_z, label=1, motif_filter=MOTIF)
    
    if X_a is None or X_z is None:
        raise ValueError("One or both input files resulted in empty data!")
        
    X = np.vstack([X_a, X_z])
    y = np.concatenate([y_a, y_z])
    
    # Generate Feature Names
    feat_names = []
    for i in range(len_a):
        feat_names.extend([f"Seq_{i+1}_A", f"Seq_{i+1}_C", f"Seq_{i+1}_G", f"Seq_{i+1}_T"])
    
    feat_names += [f"IPD_{i+1}" for i in range(len_a)] + \
                  [f"PW_{i+1}" for i in range(len_a)]
                  
    print(f"Total samples merged: {X.shape[0]} (A: {len(y_a)}, Z: {len(y_z)})")
    
    return X, y, feat_names

def sample_balanced_subset(X, y, size_per_class, seed):
    rng = np.random.RandomState(seed)
    idx_0 = np.where(y == 0)[0]
    idx_1 = np.where(y == 1)[0]
    n_0 = min(len(idx_0), size_per_class)
    n_1 = min(len(idx_1), size_per_class)
    sub_0 = rng.choice(idx_0, n_0, replace=False)
    sub_1 = rng.choice(idx_1, n_1, replace=False)
    indices = np.concatenate([sub_0, sub_1])
    rng.shuffle(indices)
    return X[indices], y[indices]

# --- PYTORCH MODEL (LOGITS OUTPUT) ---

class MLPModel(nn.Module):
    def __init__(self, input_dim, hidden_dim=128):
        super(MLPModel, self).__init__()
        # 2 Hidden Layers (128 units each)
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.bn1 = nn.BatchNorm1d(hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, hidden_dim)
        self.bn2 = nn.BatchNorm1d(hidden_dim)
        self.fc3 = nn.Linear(hidden_dim, 1)
        self.relu = nn.ReLU()
        # No sigmoid here (Logits output) for SHAP stability
        self.dropout = nn.Dropout(0.2)

    def forward(self, x):
        x = self.fc1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.dropout(x)
        
        x = self.fc2(x)
        x = self.bn2(x)
        x = self.relu(x)
        x = self.dropout(x)
        
        x = self.fc3(x) # Return logits
        return x

def train_model(model, X_train, y_train):
    model = model.to(device)
    # Use BCEWithLogitsLoss for numerical stability with logits
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=5, verbose=True)
    
    dataset = TensorDataset(torch.FloatTensor(X_train), torch.FloatTensor(y_train).unsqueeze(1))
    loader = DataLoader(dataset, batch_size=2048, shuffle=True, pin_memory=(device.type == 'cuda'))
    
    model.train()
    best_loss = float('inf')
    patience_counter = 0
    max_patience = 15 
    
    loss_history = []

    for epoch in range(200):
        epoch_loss = 0.0
        for X_batch, y_batch in loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            optimizer.zero_grad()
            outputs = model(X_batch)
            loss = criterion(outputs, y_batch)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
        
        avg_loss = epoch_loss / len(loader)
        loss_history.append(avg_loss)
        scheduler.step(avg_loss)
        
        if (epoch+1) % 10 == 0:
            print(f"Epoch {epoch+1}/200 - Loss: {avg_loss:.4f}")
            
        if avg_loss < (best_loss - 0.001):
            best_loss = avg_loss
            patience_counter = 0
        else:
            patience_counter += 1
            
        if patience_counter >= max_patience:
            print(f"Early stopping at epoch {epoch+1}")
            break
            
    return model, loss_history

# --- MAIN EXECUTION ---

print("Step 1: Loading Data...")
X, y, feat_names = load_and_merge_files(FILE_A, FILE_Z)

X_train_full, X_test, y_train_full, y_test = train_test_split(
    X, y, test_size=0.3, random_state=RANDOM_SEED, stratify=y
)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train_full)
X_test_scaled = scaler.transform(X_test)

print("\nStep 2: Training PyTorch MLP Model...")
X_sub, y_sub = sample_balanced_subset(X_train_scaled, y_train_full, 
                                      TRAIN_SAMPLE_SIZE, RANDOM_SEED)

input_dim = X_sub.shape[1]
model = MLPModel(input_dim, hidden_dim=128)
model, loss_curve = train_model(model, X_sub, y_sub)

# Save Loss Curve
df_loss = pd.DataFrame({
    'Epoch': range(1, len(loss_curve) + 1),
    'Loss': loss_curve
})
df_loss.to_csv(f'{MOTIF}_k{KMER_LENGTH}.mlp_training_loss_curve.tsv', sep='\t', index=False)

print("\nStep 3: Calculating SHAP Values (Using GradientExplainer)...")

# GradientExplainer is generally more robust for PyTorch than DeepExplainer
# We must ensure the model is in eval mode
model.eval()

# Prepare Background Data
background_size = 500
background_idx = np.random.choice(X_sub.shape[0], background_size, replace=False)
background_data = torch.FloatTensor(X_sub[background_idx]).to(device)

# Prepare Test Data
test_size = 500 
X_test_sample_np = X_test_scaled[:test_size]
X_test_sample = torch.FloatTensor(X_test_sample_np).to(device)
X_test_sample.requires_grad = True # Required for GradientExplainer

try:
    # Try GradientExplainer first
    explainer = shap.GradientExplainer(model, background_data)
    shap_values_raw = explainer.shap_values(X_test_sample)
except Exception as e:
    print(f"GradientExplainer failed: {e}")
    print("Falling back to DeepExplainer...")
    try:
        explainer = shap.DeepExplainer(model, background_data)
        shap_values_raw = explainer.shap_values(X_test_sample)
    except Exception as e2:
        print(f"DeepExplainer failed: {e2}")
        print("CRITICAL FAIL: Neither Gradient nor Deep explainer worked. Exiting.")
        sys.exit(1)

# Robustly handle list or array output
shap_values_final = None

if isinstance(shap_values_raw, list):
    # Depending on SHAP version and model output (single unit), list might contain one array
    # If binary classification (single output), shap usually returns single array or list of one array
    shap_values_final = np.array(shap_values_raw[0]) 
else:
    shap_values_final = np.array(shap_values_raw)

# Fix: Ensure shap_values is a numpy array and correct dimensions
if torch.is_tensor(shap_values_final):
    shap_values_final = shap_values_final.cpu().detach().numpy()

# Squeeze if necessary (remove dimension of size 1)
shap_values_final = np.squeeze(shap_values_final)

# Ensure shape matches (N_samples, N_features)
if shap_values_final.ndim != 2:
    print(f"Warning: Unexpected SHAP values shape: {shap_values_final.shape}. Reshaping...")
    # Try to reshape to (test_size, num_features)
    try:
        shap_values_final = shap_values_final.reshape(test_size, -1)
    except:
        print("Error: Could not reshape SHAP values.")
        sys.exit(1)

print(f"Final SHAP matrix shape: {shap_values_final.shape}")

print("Saving SHAP results...")
# Summary DataFrame
mean_abs_shap = np.mean(np.abs(shap_values_final), axis=0)
df_shap_summary = pd.DataFrame({
    'Feature': feat_names,
    'Mean_Abs_SHAP': mean_abs_shap
}).sort_values(by='Mean_Abs_SHAP', ascending=False)
df_shap_summary.to_csv(f'{MOTIF}_k{KMER_LENGTH}.shap_feature_importance_summary.tsv', sep='\t', index=False)

# Raw Values
df_shap_raw = pd.DataFrame(shap_values_final, columns=feat_names)
df_features_raw = pd.DataFrame(X_test_sample_np, columns=[f"{f}_Value" for f in feat_names])
df_shap_full = pd.concat([df_features_raw, df_shap_raw], axis=1)
df_shap_full.to_csv(f'{MOTIF}_k{KMER_LENGTH}.shap_values_subset_for_beeswarm.tsv', sep='\t', index=False)

# Plotting
plt.figure(figsize=(10, 8))
shap.summary_plot(shap_values_final, X_test_sample_np, feature_names=feat_names, 
                 show=False, plot_type="bar", color='#44aa99')
plt.title(f"Feature Importance ({KMER_LENGTH}-mer)", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{MOTIF}_k{KMER_LENGTH}.shap_feature_importance.pdf', format='pdf')
plt.close()

plt.figure(figsize=(10, 8))
shap.summary_plot(shap_values_final, X_test_sample_np, feature_names=feat_names, 
                 show=False, cmap='viridis') 
plt.title(f"SHAP Summary ({MOTIF}, {KMER_LENGTH}-mer)", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{MOTIF}_k{KMER_LENGTH}.shap_beeswarm.pdf', format='pdf')
plt.close()

print("\nAnalysis Complete.")
