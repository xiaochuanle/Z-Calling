#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import pysam
import re

def reverse_complement(seq):
    # Function to compute the reverse complement of a DNA sequence
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in reversed(seq)])

def get_context_kmer(read, pos, base, k):
    # Extracts the k-mer context for the target base
    half_k = (k - 1) // 2
    if base == 'A':
        return read.query_sequence[(pos-half_k):(pos+half_k+1)]
    elif base == 'T':
        kmer = read.query_sequence[(pos-half_k):(pos+half_k+1)]
        return reverse_complement(kmer)

def get_kinetics(read, pos, base, k):
    # Retrieves the kinetics information (fi, fp, ri, rp) for the context bases
    half_k = (k - 1) // 2
    if base == 'A':
        ipd = read.get_tag('ri')[(len(read.query_sequence)-1-pos-half_k):(len(read.query_sequence)-1-pos+half_k+1)]
        pw = read.get_tag('rp')[(len(read.query_sequence)-1-pos-half_k):(len(read.query_sequence)-1-pos+half_k+1)]
    elif base == 'T':
        ipd = read.get_tag('fi')[(pos-half_k):(pos+half_k+1)]
        pw = read.get_tag('fp')[(pos-half_k):(pos+half_k+1)]
    return ipd, pw

def process_bam(bam_stream, k, regex=None):
    results = []
    pattern = re.compile(regex) if regex else None
    half_k = (k - 1) // 2
    pattern_len = len(regex) if regex else 0
    half_pattern = (pattern_len - 1) // 2

    with pysam.AlignmentFile(bam_stream, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            # Skip reads without sequence
            if read.query_sequence is None:
                continue

            # Check for the presence and length of tags before processing
            has_fi_fp = len(read.get_tag('fi')) == len(read.query_sequence) and len(read.get_tag('fp')) == len(read.query_sequence)
            has_ri_rp = len(read.get_tag('ri')) == len(read.query_sequence) and len(read.get_tag('rp')) == len(read.query_sequence)
            
            if not has_fi_fp or not has_ri_rp:
                continue

            for pos, base in enumerate(read.query_sequence):
                if base not in ['A', 'T']:
                    continue
                if pos-half_k < 0 or pos+half_k >= len(read.query_sequence):
                    continue
                if pos < 50 or pos > len(read.query_sequence)-50:
                    continue
                context_kmer = get_context_kmer(read, pos, base, k)
                center_start = half_k - half_pattern
                center_end = half_k + half_pattern + 1
                center_kmer = context_kmer[center_start:center_end]
                if pattern and not pattern.match(center_kmer):
                    continue
                if context_kmer.count('A') - center_kmer.count('A') > 0:
                    continue

                kinetics = get_kinetics(read, pos, base, k)
                ipd = kinetics[0].tolist()
                pw = kinetics[1].tolist()
                results.append([pos, base, context_kmer, ipd, pw])

    columns = ['Position', 'Base', 'ContextKmer', 'ContextIPD', 'ContextPW']
    return pd.DataFrame(results, columns=columns)

def main():
    parser = argparse.ArgumentParser(description="Process a BAM file and extract sequencing information.")
    parser.add_argument('bam_file', type=str, help="Path to the BAM file.")
    parser.add_argument('-o', '--output', type=str, default=None, help="Path to output CSV file. If not specified, prints to stdout.")
    parser.add_argument('-k', '--kmer_length', type=int, default=21, help="Length of target base context to be output.")
    parser.add_argument('-r', '--regex', type=str, default=None, help="Regex pattern to filter kmers based on context sequence.")

    args = parser.parse_args()
    k_length = args.kmer_length

    df = process_bam(args.bam_file, k_length, args.regex)

    if args.output:
        df.to_csv(args.output, sep="\t", index=False)
        print(f"Results saved to {args.output}")
    else:
        print(df)

if __name__ == "__main__":
    main()
