#!/usr/bin/env python3
import argparse
import pysam
import numpy as np
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Filter BAM/SAM reads based on median ratios of tags.')
    parser.add_argument('-b', '--bam', type=str, required=True, help='Path to input BAM/SAM file.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to output BAM file.')
    return parser.parse_args()

def compute_median_ratio(read):
    # Extract the tags
    try:
        fi = read.get_tag('fi')
        ri = read.get_tag('ri')
        fp = read.get_tag('fp')
        rp = read.get_tag('rp')
    except KeyError:
        # If any tag is missing, return None to indicate the read should be skipped
        return None

    # Ensure tag values are lists or arrays
    for tag_name, tag_values in [('fi', fi), ('ri', ri), ('fp', fp), ('rp', rp)]:
        if not isinstance(tag_values, (list, tuple, np.ndarray)):
            tag_values = [tag_values]
        # Convert to numpy array for computation
        tag_values_array = np.array(tag_values)
        # Check for empty arrays
        if tag_values_array.size == 0:
            return None  # Skip reads with empty tag values

        # Store back the arrays
        if tag_name == 'fi':
            fi_array = tag_values_array
        elif tag_name == 'ri':
            ri_array = tag_values_array
        elif tag_name == 'fp':
            fp_array = tag_values_array
        elif tag_name == 'rp':
            rp_array = tag_values_array

    # Compute median for each tag
    medians = {
        'fi': np.median(fi_array),
        'ri': np.median(ri_array),
        'fp': np.median(fp_array),
        'rp': np.median(rp_array)
    }

    # Compute median_ratio
    numerator = medians['fp'] + medians['rp']
    denominator = medians['fi'] + medians['ri']
    if denominator == 0:
        return None  # Avoid division by zero; skip this read
    median_ratio = numerator / denominator

    return median_ratio

def check_initial_reads(bam_file):
    # Check the first 50 reads for required tags
    bam_in = pysam.AlignmentFile(bam_file, 'rb', check_sq=False)
    count = 0
    for read in bam_in.fetch(until_eof=True):
        try:
            read.get_tag('fi')
            read.get_tag('ri')
            read.get_tag('fp')
            read.get_tag('rp')
        except KeyError:
            pass
        else:
            count += 1
        if count > 0:
            break
        if bam_in.tell() >= 50:
            break

    bam_in.close()

    if count == 0:
        print("BAM does not contain kinetic information required for Z-Calling. Exit!")
        sys.exit(1)

def main():
    args = parse_args()

    # Check initial reads for required tags
    check_initial_reads(args.bam)

    # Open input BAM/SAM file
    bam_file = args.bam
    if bam_file.endswith('.bam') or bam_file.endswith('.BAM'):
        bam_in = pysam.AlignmentFile(bam_file, 'rb', check_sq=False)
    else:
        bam_in = pysam.AlignmentFile(bam_file, 'r', check_sq=False)

    # First pass: Compute median_ratios for all reads
    median_ratios = []
    for read in bam_in.fetch(until_eof=True):
        # Skip reads without the necessary tags
        try:
            fn = read.get_tag('fn')
            rn = read.get_tag('rn')
        except KeyError:
            continue

        if fn < 3 or rn < 3:
            continue  # Skip reads with fn < 3 or rn < 3

        median_ratio = compute_median_ratio(read)
        if median_ratio is not None:
            median_ratios.append(median_ratio)
        # Reads with median_ratio == None are skipped automatically

    bam_in.close()

    if not median_ratios:
        print("No reads with valid median_ratios found. Exiting.")
        sys.exit(1)

    # Compute thresholds
    median_ratios_array = np.array(median_ratios)
    percentile_99 = np.percentile(median_ratios_array, 99)

    # Second pass: Filter reads and write to output BAM file
    bam_in = pysam.AlignmentFile(bam_file, 'rb', check_sq=False)
    bam_out = pysam.AlignmentFile(args.output, 'wb', template=bam_in)

    for read in bam_in.fetch(until_eof=True):
        # Skip reads without the necessary tags
        try:
            fn = read.get_tag('fn')
            rn = read.get_tag('rn')
        except KeyError:
            continue

        if fn < 3 or rn < 3:
            continue  # Skip reads with fn < 3 or rn < 3

        median_ratio = compute_median_ratio(read)
        if median_ratio is None:
            continue  # Skip reads missing any of the four tags
        if median_ratio < 0.3:
            continue  # Exclude reads with median_ratio < 0.3
        if median_ratio >= percentile_99:
            continue  # Exclude top 1% reads with largest median_ratios
        # Write read to output BAM file
        bam_out.write(read)

    bam_in.close()
    bam_out.close()

    print(f"Filtering complete. Output written to {args.output}")

if __name__ == '__main__':
    main()
