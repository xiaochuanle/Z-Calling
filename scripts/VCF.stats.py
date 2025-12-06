#!/usr/bin/env python3
import sys
from collections import defaultdict

rev_comp = {"A->T": "T->A", "A->G": "T->C", "A->C": "T->G", "C->A": "G->T", "C->T": "G->A", "C->G": "G->C"}

def reverse_complement(substitution):
    return rev_comp.get(substitution, substitution)

def process_vcf(vcf_file):
    substitution_counts = defaultdict(int)
    with open(vcf_file, 'r') as file:
        for line in file:

            # Extract REF and ALT columns
            parts = line.strip().split('\t')
            ref = parts[3]
            alt_list = parts[4].split(',')
            if ref == "N":
                continue

            # Ignore non-sense alleles
            alt_list = [alt for alt in alt_list if alt != '<*>']

            for alt in alt_list:
                if len(ref) == 1 and len(alt) == 1:
                    # Process single base substitutions
                    substitution = f"{ref}->{alt}"
                    if substitution in rev_comp:
                        rev_substitution = reverse_complement(substitution)
                        substitution_counts[rev_substitution] += 1
                    else:
                        substitution_counts[substitution] += 1
                else:
                    # Process INDELs
                    if len(alt) > len(ref):
                        substitution_counts['INS'] += 1
                    elif len(alt) < len(ref):
                        substitution_counts['DEL'] += 1
                    else:
                        continue  # Shouldn't happen, but just in case

    return substitution_counts

def save_to_tsv(substitution_counts, output_file):
    with open(output_file, 'w') as out_file:
        out_file.write("Substitution_Type\tCount\n")
        for substitution, count in substitution_counts.items():
            out_file.write(f"{substitution}\t{count}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_vcf> <output_tsv>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_tsv = sys.argv[2]

    substitution_counts = process_vcf(input_vcf)
    sorted_dict = {key: substitution_counts[key] for key in sorted(substitution_counts.keys())}
    save_to_tsv(sorted_dict, output_tsv)
