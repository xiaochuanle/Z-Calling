#!/usr/bin/env python3
import argparse
import pysam
import numpy as np
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Process BAM/SAM file and compute statistics for fi, ri, fp, rp tags.')
    parser.add_argument('-b', '--bam', type=str, help='Path to BAM/SAM file. Use "-" or omit for SAM input from stdin.')
    parser.add_argument('-s', '--status', type=str, required=True, help='Path to status TSV file.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output TSV file name.')
    return parser.parse_args()

def read_status_file(status_file):
    status_dict = {}
    with open(status_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 2:
                continue
            read_name, status = parts
            status_dict[read_name] = status
    return status_dict

def process_read(read, status_dict):
    read_name = read.query_name

    # Extract the tags
    try:
        fi = read.get_tag('fi')
        ri = read.get_tag('ri')
        fp = read.get_tag('fp')
        rp = read.get_tag('rp')
    except KeyError:
        # If any tag is missing, skip this read
        return None

    # Compute statistics for each tag
    result = [read_name]
    for tag_values in [fi, ri, fp, rp]:
        if not isinstance(tag_values, (list, tuple, np.ndarray)):
            tag_values = [tag_values]
        tag_values = np.array(tag_values)

        # Handle empty tag_values
        if tag_values.size == 0:
            median = 'NA'
            q25 = 'NA'
            q75 = 'NA'
            q90 = 'NA'
        else:
            median = np.median(tag_values)
            q25 = np.percentile(tag_values, 25)
            q75 = np.percentile(tag_values, 75)
            q90 = np.percentile(tag_values, 90)
        result.extend([median, q25, q75, q90])

    # Append the status
    status = status_dict.get(read_name, 'NA')
    result.append(status)

    return result

def main():
    args = parse_args()
    status_dict = read_status_file(args.status)

    # Open BAM or SAM file with check_sq=False
    if args.bam and args.bam != '-':
        bam_file = args.bam
        if bam_file.endswith('.bam') or bam_file.endswith('.BAM'):
            bam = pysam.AlignmentFile(bam_file, 'rb', check_sq=False)
        else:
            bam = pysam.AlignmentFile(bam_file, 'r', check_sq=False)
    else:
        bam = pysam.AlignmentFile('-', 'r', check_sq=False)

    # Prepare the output file
    output_file = open(args.output, 'w')

    # Process reads one by one
    for read in bam.fetch(until_eof=True):
        result = process_read(read, status_dict)
        if result and len(result) == 18:
            output_file.write('\t'.join(map(str, result)) + '\n')

    bam.close()
    output_file.close()

if __name__ == '__main__':
    main()

