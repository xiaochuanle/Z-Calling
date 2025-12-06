#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import os
import joblib

# --- ARGUMENT PARSING ---
parser = argparse.ArgumentParser(description="Benchmark classification models on kinetic data (CUDA enabled).")
parser.add_argument("file_a", help="Input file for Canonical A (TSV)")
parser.add_argument("file_z", help="Input file for Z-base (TSV)")
parser.add_argument("motif", help="Motif to filter by (e.g., TTATT)")
parser.add_argument("-k", "--kmer_length", type=int, default=21, 
                    help="Length of center k-mer to use (odd number, <= 21). Default: 21")
parser.add_argument("-s", "--sample_size", type=int, default=100000,
                    help="Number of samples per class for training. Default: 100000")
parser.add_argument("-b", "--batch_size", type=int, default=1024,
                    help="Batch size for training. Default: 1024")

args = parser.parse_args()

FILE_A = args.file_a
FILE_Z = args.file_z
MOTIF = args.motif
KMER_LENGTH = args.kmer_length
TRAIN_SAMPLE_SIZE = args.sample_size
BATCH_SIZE = args.batch_size

# Input sequence length (assumed from your file description)
INPUT_SEQ_LEN = 21 

if KMER_LENGTH % 2 == 0:
    print(f"Error: K-mer length must be odd. Got {KMER_LENGTH}.")
    sys.exit(1)
if KMER_LENGTH > INPUT_SEQ_LEN:
    print(f"Error: K-mer length ({KMER_LENGTH}) cannot be larger than input length ({INPUT_SEQ_LEN}).")
    sys.exit(1)

# Calculate slicing indices for center k-mer
SLICE_START = (INPUT_SEQ_LEN - KMER_LENGTH) // 2
SLICE_END = SLICE_START + KMER_LENGTH

print(f"Configuration: Motif={MOTIF}, K-mer Length={KMER_LENGTH}, Train Size={TRAIN_SAMPLE_SIZE}, Batch Size={BATCH_SIZE}")
print(f"  - Slicing input arrays from index {SLICE_START} to {SLICE_END}")

# --- DEVICE CONFIGURATION (CUDA SUPPORT) ---
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"  - PyTorch computing device: {device}")
# -------------------------------------------

# Constants
RANDOM_SEED = 42
# TRAIN_SAMPLE_SIZE set via args
# BATCH_SIZE set via args
MAX_EPOCHS = 1000  # Increased as requested
PATIENCE = 15      # Increased to avoid premature stopping
MIN_DELTA = 0.0001 # Refined delta
MODEL_DIR = f'{MOTIF}_benchmark_models'

if not os.path.exists(MODEL_DIR):
    os.makedirs(MODEL_DIR)

# --- DATA LOADING & PREPROCESSING ---
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
    # Slice the string first
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
        return None, None, None

    if motif_filter and 'ContextKmer' in df.columns:
        original_count = len(df)
        df = df[df['ContextKmer'].str.contains(motif_filter)]
        print(f"  - Filtered {original_count} -> {len(df)} rows with motif '{motif_filter}'")
    
    if df.empty: return None, None, None

    ipd_data = np.stack(df['ContextIPD'].apply(parse_kinetic_string_and_slice).values)
    pw_data = np.stack(df['ContextPW'].apply(parse_kinetic_string_and_slice).values)
    
    if 'ContextKmer' not in df.columns:
        raise ValueError(f"File {filepath} missing 'ContextKmer' column required for sequence features.")
    seq_data = np.stack(df['ContextKmer'].apply(one_hot_encode_sequence_slice).values)
    
    # Combine: Sequence + IPD + PW
    X_part = np.hstack([seq_data, ipd_data, pw_data])
    y_part = np.full(X_part.shape[0], label)
    
    return X_part, y_part, KMER_LENGTH 

def load_and_merge_files(file_a, file_z):
    X_a, y_a, len_a = process_single_file(file_a, 0, MOTIF)
    X_z, y_z, _ = process_single_file(file_z, 1, MOTIF)
    
    if X_a is None or X_z is None: raise ValueError("Empty input files")
        
    X = np.vstack([X_a, X_z])
    y = np.concatenate([y_a, y_z])
    
    input_dim_total = X.shape[1] 
    print(f"Total samples: {X.shape[0]}. Input Feature Dimension: {input_dim_total}")
    
    return X, y, input_dim_total

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

# --- MODEL DEFINITIONS (PyTorch) ---

class MLPModel(nn.Module):
    def __init__(self, input_dim, hidden_dim=128): # Default changed to 128
        super(MLPModel, self).__init__()
        # 2 Hidden Layers (128 units each)
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.bn1 = nn.BatchNorm1d(hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, hidden_dim)
        self.bn2 = nn.BatchNorm1d(hidden_dim)
        self.fc3 = nn.Linear(hidden_dim, 1)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
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
        
        x = self.fc3(x)
        x = self.sigmoid(x)
        return x

class CNNModel(nn.Module):
    def __init__(self, total_features):
        super(CNNModel, self).__init__()
        self.kmer_len = KMER_LENGTH
        self.num_channels = 6 
        
        self.conv1 = nn.Conv1d(in_channels=self.num_channels, out_channels=32, kernel_size=3, padding=1)
        self.relu = nn.ReLU()
        pool_k = 2 if self.kmer_len >= 2 else 1
        self.pool = nn.MaxPool1d(kernel_size=pool_k)
        
        self.conv2 = nn.Conv1d(32, 64, kernel_size=3, padding=1)
        self.flatten = nn.Flatten()
        
        l_out = self.kmer_len // pool_k // pool_k
        flat_size = l_out * 64
        
        self.fc1 = nn.Linear(flat_size, 64)
        self.fc2 = nn.Linear(64, 1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        batch_size = x.size(0)
        seq = x[:, :self.kmer_len*4].view(batch_size, self.kmer_len, 4)
        ipd = x[:, self.kmer_len*4:self.kmer_len*5].view(batch_size, self.kmer_len, 1)
        pw = x[:, self.kmer_len*5:].view(batch_size, self.kmer_len, 1)
        
        x_reshaped = torch.cat([seq, ipd, pw], dim=2).permute(0, 2, 1)
        
        x = self.pool(self.relu(self.conv1(x_reshaped)))
        x = self.pool(self.relu(self.conv2(x)))
        x = self.flatten(x)
        x = self.relu(self.fc1(x))
        x = self.sigmoid(self.fc2(x))
        return x

class BiLSTMModel(nn.Module):
    def __init__(self, total_features, hidden_dim=64):
        super(BiLSTMModel, self).__init__()
        self.kmer_len = KMER_LENGTH
        self.lstm = nn.LSTM(input_size=6, hidden_size=hidden_dim, batch_first=True, bidirectional=True)
        self.fc = nn.Linear(hidden_dim * 2, 1) 
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        batch_size = x.size(0)
        seq = x[:, :self.kmer_len*4].view(batch_size, self.kmer_len, 4)
        ipd = x[:, self.kmer_len*4:self.kmer_len*5].view(batch_size, self.kmer_len, 1)
        pw = x[:, self.kmer_len*5:].view(batch_size, self.kmer_len, 1)
        
        x_reshaped = torch.cat([seq, ipd, pw], dim=2)
        
        lstm_out, _ = self.lstm(x_reshaped)
        last_hidden = lstm_out[:, -1, :] 
        out = self.sigmoid(self.fc(last_hidden))
        return out

class GCNModel(nn.Module):
    def __init__(self, total_features):
        super(GCNModel, self).__init__()
        self.kmer_len = KMER_LENGTH
        self.gcn1 = nn.Conv1d(6, 64, kernel_size=3, padding=1) 
        self.gcn2 = nn.Conv1d(64, 64, kernel_size=3, padding=1)
        self.fc = nn.Linear(self.kmer_len * 64, 1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        batch_size = x.size(0)
        seq = x[:, :self.kmer_len*4].view(batch_size, self.kmer_len, 4)
        ipd = x[:, self.kmer_len*4:self.kmer_len*5].view(batch_size, self.kmer_len, 1)
        pw = x[:, self.kmer_len*5:].view(batch_size, self.kmer_len, 1)
        x_reshaped = torch.cat([seq, ipd, pw], dim=2).permute(0, 2, 1)
        
        x = torch.relu(self.gcn1(x_reshaped))
        x = torch.relu(self.gcn2(x))
        x = x.view(x.size(0), -1)
        x = self.sigmoid(self.fc(x))
        return x

def train_pytorch_model(model, X_train, y_train, max_epochs=MAX_EPOCHS, patience=PATIENCE, min_delta=MIN_DELTA):
    model = model.to(device)
    criterion = nn.BCELoss()
    # LOWER LEARNING RATE for stability
    optimizer = optim.Adam(model.parameters(), lr=0.0001) 
    
    # NEW: Learning Rate Scheduler (Reduce LR if loss plateaus)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=3, verbose=True)

    dataset = TensorDataset(torch.FloatTensor(X_train), torch.FloatTensor(y_train).unsqueeze(1))
    loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True, pin_memory=(device.type == 'cuda'))
    
    best_loss = float('inf')
    epochs_no_improve = 0
    stopped_early = False
    
    model.train()
    for epoch in range(max_epochs):
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
        
        # Step the scheduler
        scheduler.step(avg_loss)
        
        if (epoch + 1) % 10 == 0:
            print(f"    Epoch {epoch+1}/{max_epochs} - Loss: {avg_loss:.4f}")

        if avg_loss < (best_loss - min_delta):
            best_loss = avg_loss
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1
            
        if epochs_no_improve >= patience:
            print(f"    - Early stopping at epoch {epoch+1} (Best Loss: {best_loss:.4f})")
            stopped_early = True
            break
            
    if not stopped_early:
        print(f"    - Completed {max_epochs} epochs. Loss: {avg_loss:.4f}.")
        
    return model

def evaluate_pytorch_model(model, X_test):
    model.eval()
    # Create a DataLoader for the test set to batch evaluation
    dataset = TensorDataset(torch.FloatTensor(X_test))
    # Use a reasonable batch size for inference (can be larger than training batch size)
    loader = DataLoader(dataset, batch_size=BATCH_SIZE * 2, shuffle=False)
    
    all_preds = []
    
    with torch.no_grad():
        for X_batch in loader:
            # Move batch to GPU
            X_batch = X_batch[0].to(device)
            # Predict
            preds = model(X_batch).cpu().squeeze().numpy()
            # Handle case where batch size is 1 (squeeze might remove too many dims)
            if preds.ndim == 0:
                preds = np.array([preds])
            all_preds.append(preds)
            
    return np.concatenate(all_preds)

# --- MAIN BENCHMARK FUNCTION ---

def run_benchmark():
    print(f"Starting benchmark for Motif: {MOTIF}, K-mer: {KMER_LENGTH}")
    X, y, input_dim = load_and_merge_files(FILE_A, FILE_Z)
    
    X_train_full, X_test, y_train_full, y_test = train_test_split(
        X, y, test_size=0.3, random_state=RANDOM_SEED, stratify=y
    )
    
    scaler = StandardScaler()
    X_train_full_scaled = scaler.fit_transform(X_train_full)
    X_test_scaled = scaler.transform(X_test)
    
    models = ['RandomForest', 'MLP', 'CNN', 'BiLSTM', 'GCN']
    seeds = [42, 123, 999]
    results = []

    print(f"Benchmarking {len(models)} models x {len(seeds)} replicates...")

    for model_name in models:
        for seed in seeds:
            model_path = os.path.join(MODEL_DIR, f'{model_name}_seed{seed}.pkl')
            
            if os.path.exists(model_path):
                print(f"  Loading pre-trained {model_name} (Seed {seed})...")
                if model_name in ['RandomForest', 'MLP']:
                    clf = joblib.load(model_path)
                    preds = clf.predict_proba(X_test_scaled)[:, 1]
                else:
                    if model_name == 'MLP': model = MLPModel(input_dim)
                    elif model_name == 'CNN': model = CNNModel(input_dim)
                    elif model_name == 'BiLSTM': model = BiLSTMModel(input_dim)
                    elif model_name == 'GCN': model = GCNModel(input_dim)
                    
                    model.load_state_dict(torch.load(model_path, map_location=device))
                    model = model.to(device)
                    preds = evaluate_pytorch_model(model, X_test_scaled)
            else:
                print(f"  Training {model_name} (Seed {seed})...")
                X_sub, y_sub = sample_balanced_subset(X_train_full_scaled, y_train_full, TRAIN_SAMPLE_SIZE, seed)
                preds = None
                
                if model_name == 'RandomForest':
                    clf = RandomForestClassifier(n_estimators=100, random_state=seed, n_jobs=-1)
                    clf.fit(X_sub, y_sub)
                    joblib.dump(clf, model_path)
                    preds = clf.predict_proba(X_test_scaled)[:, 1]
                    
                elif model_name == 'MLP':
                    # Use PyTorch MLP
                    model = MLPModel(input_dim)
                    train_pytorch_model(model, X_sub, y_sub)
                    torch.save(model.state_dict(), model_path)
                    preds = evaluate_pytorch_model(model, X_test_scaled)
                    
                elif model_name == 'CNN':
                    model = CNNModel(input_dim)
                    train_pytorch_model(model, X_sub, y_sub)
                    torch.save(model.state_dict(), model_path)
                    preds = evaluate_pytorch_model(model, X_test_scaled)
                    
                elif model_name == 'BiLSTM':
                    model = BiLSTMModel(input_dim) 
                    train_pytorch_model(model, X_sub, y_sub)
                    torch.save(model.state_dict(), model_path)
                    preds = evaluate_pytorch_model(model, X_test_scaled)
                    
                elif model_name == 'GCN':
                    model = GCNModel(input_dim)
                    train_pytorch_model(model, X_sub, y_sub)
                    torch.save(model.state_dict(), model_path)
                    preds = evaluate_pytorch_model(model, X_test_scaled)
            
            score = roc_auc_score(y_test, preds)
            results.append({'Model': model_name, 'Replicate': seed, 'AUC': score})

    output_tsv = f"{MOTIF}_k{KMER_LENGTH}_model_benchmark_results.tsv"
    df_res = pd.DataFrame(results)
    df_res.to_csv(output_tsv, sep='\t', index=False)
    print(f"Results saved to {output_tsv}")
    
    output_pdf = f"{MOTIF}_k{KMER_LENGTH}_model_benchmark_boxplot.pdf"
    plt.figure(figsize=(10, 6))
    
    sns.boxplot(x='Model', y='AUC', hue='Model', data=df_res, palette="viridis", legend=False)
    sns.swarmplot(x='Model', y='AUC', data=df_res, color=".25")
    
    plt.title(f'Performance Comparison ({MOTIF}, {KMER_LENGTH}-mer)', fontsize=14)
    plt.ylabel('AUC Score')
    plt.ylim(0.8, 1.0) 
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.savefig(output_pdf, format='pdf')
    print(f"Plot saved to {output_pdf}")

if __name__ == "__main__":
    run_benchmark()
