#!/usr/bin/env python3

import pandas as pd
import sys
import argparse

def main():
    # Setting up argument parser
    parser = argparse.ArgumentParser(description='Process SVM classifier output results.')
    parser.add_argument('--SVM', required=True, help='Path to the SVM output file (mandatory)')
    parser.add_argument('--read_ref', required=True, help='Path to the read-to-RefContig mapping file (mandatory)')
    parser.add_argument('--ref_species', required=False, help='Path to the RefContig-to-Species mapping file (optional)')
    parser.add_argument('--output', required=True, help='Path to the output TSV file (mandatory)')
    args = parser.parse_args()

    # Reading the input files
    SVM = pd.read_table(args.SVM, names=['read', 'tag'])
    READ = pd.read_table(args.read_ref, names=['read', 'RefContig'])

    # Filter out reads with tag '*'
    SVM = SVM.loc[SVM['tag'] != '*']
    SVM = SVM.merge(READ, on='read').copy()

    # If ref_species is given, merge with TAXON to include 'Species' information
    if args.ref_species:
        TAXON = pd.read_table(args.ref_species, names=['RefContig', 'Species'])
        SVM = SVM.merge(TAXON, on='RefContig').copy()[['read', 'tag', 'Species']]
        group_col = 'Species'
    else:
        # Otherwise, use 'RefContig' as the grouping column
        SVM = SVM[['read', 'tag', 'RefContig']]
        group_col = 'RefContig'

    # Calculate read counts grouped by the specified column (either 'Species' or 'RefContig')
    SVM_SpeciesReadCount = SVM[['read', group_col]].groupby(group_col).count().copy()
    SVM_SpeciesZreadCount = SVM.loc[SVM['tag'] == '+'][['read', group_col]].groupby(group_col).count().copy()

    # Calculate likelihood ratio (LR+)
    SVM_LR = SVM_SpeciesZreadCount / SVM_SpeciesReadCount / 0.002

    # Reset indices and merge read counts
    SVM_SpeciesReadCount = SVM_SpeciesReadCount.reset_index()
    SVM_SpeciesZreadCount = SVM_SpeciesZreadCount.reset_index()
    MERGE = SVM_SpeciesReadCount.merge(SVM_SpeciesZreadCount, on=group_col, how='left')

    # Fill missing values with 0
    MERGE.fillna(0, inplace=True)

    # Rename columns for clarity and calculate LR+
    MERGE['ClassifiedReadCt'] = MERGE.read_x.astype(int)
    MERGE['ClassifiedZReadCt'] = MERGE.read_y.astype(int)
    MERGE['LR+'] = MERGE.ClassifiedZReadCt / MERGE.ClassifiedReadCt / 0.002

    # Drop unnecessary columns and save the output to a TSV file
    MERGE.drop(['read_x', 'read_y'], axis=1).to_csv(args.output, sep="\t", index=False, header=True)

if __name__ == "__main__":
    main()
