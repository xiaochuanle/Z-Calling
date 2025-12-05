#!/bin/bash

# 1. Setup Environment
# Add the build directory to the PATH so we can run z-calling-base, z-bam2txt, etc.
# Note: Ensure you have activated the conda environment before running this script:
# conda activate Z-Calling
export PATH=/path/to/Z-Calling/build:$PATH

# 2. Filter BAM Files
# Pre-process the raw BAM files to filter out low-quality reads or artifacts.
# Processing the dZ-DNA (modified) sample:
python3 /path/to/Z-Calling/py/filter_bam.py -b Ecoli.dZ-DNA.100.bam -o Ecoli.dZ-DNA.100.filt.bam
# Processing the native DNA (control) sample:
python3 /path/to/Z-Calling/py/filter_bam.py -b Ecoli.nativeDNA.100.bam -o Ecoli.nativeDNA.100.filt.bam

# 3. Align Reads to Reference (using pbmm2)
# Align the filtered reads to the E. coli reference genome. 
# --preset CCS: Optimized for PacBio HiFi/CCS reads.
# --sort: Automatically sorts the output BAM (required for downstream steps).
pbmm2 align Ref/ECOLI.reference.fasta Ecoli.dZ-DNA.100.filt.bam Ecoli.dZ-DNA.100.filt.srt.bam --preset CCS --sort
pbmm2 align Ref/ECOLI.reference.fasta Ecoli.nativeDNA.100.filt.bam Ecoli.nativeDNA.100.filt.srt.bam --preset CCS --sort

# 4. Read-Level Classification (running k21-full-ZA MLP model and SVM classification)
# Perform base-level modification calling using the K21 model trained for read classification.
# -k 21: Use k-mer size 21, must matching the applied model.
# -t 16: Use 16 threads for parallel processing.
# Format: z-calling-base [options] <input_bam> <model_path> <output_bam>

# Run on Native DNA (Control):
z-calling-base -k 21 -t 16 Ecoli.nativeDNA.100.filt.srt.bam /path/to/Z-Calling/model/k21-full-ZA/scripted_m21.pth Ecoli.nativeDNA.100.filt.srt.k21Z.bam
# Run on dZ-DNA (Sample):
z-calling-base -k 21 -t 16 Ecoli.dZ-DNA.100.filt.srt.bam  /path/to/Z-Calling/model/k21-full-ZA/scripted_m21.pth Ecoli.dZ-DNA.100.filt.srt.k21Z.bam

# Predict the modification status of entire reads using a pre-trained SVM classifier.
# Format: z-calling-read predict <input_bam> <output_tsv> <model_path> <max_reads> <min_len>
z-calling-read predict Ecoli.dZ-DNA.100.filt.srt.k21Z.bam Ecoli.dZ-DNA.100.filt.srt.k21Z.SVM.tsv /path/to/Z-Calling/model/svm/k21ReadClassifier 500 100
z-calling-read predict Ecoli.nativeDNA.100.filt.srt.k21Z.bam Ecoli.nativeDNA.100.filt.srt.k21Z.SVM.tsv /path/to/Z-Calling/model/svm/k21ReadClassifier 500 100

# 5. A/Z Base Calling (K11 Model)
# Demonstrate calling with a different k-mer size (K11) on the dZ-DNA sample.
z-calling-base -k 11 -t 16 Ecoli.dZ-DNA.100.filt.srt.bam  /path/to/Z-Calling/model/k11-mixed-AZ/scripted_m11.pth Ecoli.dZ-DNA.100.filt.srt.k11AZ.bam

# 6. Sequence Extraction
# Extract FASTA sequences from the K11 processed BAM file.
# Z-Calling distinguishes between T-A and T-Z pairs, generating output in a six-alphabet FASTA format, where 'O' represents a T paired with a Z.
z-seq Ecoli.dZ-DNA.100.filt.srt.k11AZ.bam Ecoli.dZ-DNA.100.filt.k11AZ.fasta

# 8. Frequency Analysis (Site-Level Aggregation)
# Step A: Convert the annotated BAM to a TSV text format containing modification probabilities.
# -zp: Use the 'ZP' tag for probabilities.
z-bam2txt -zp -t 16 Ecoli.dZ-DNA.100.filt.srt.k11AZ.bam - Ref/ECOLI.reference.fasta > Ecoli.dZ-DNA.100.filt.k11AZ.tsv

# Step B: Calculate aggregated Z-base frequencies across reference coordinates.
# This generates the final frequency report (chr, coordinate, num_records_above_threshold, num_records_below_threshold, Z ratio).
zfreq -i Ecoli.dZ-DNA.100.filt.k11AZ.tsv -o Ecoli.dZ-DNA.100.filt.k11AZ.freq.tsv
