#!/bin/bash

# Define usage function to provide help information
usage() {
  echo "Usage: Z-Calling <subcommand> [options]"
  echo "Subcommands:"
  echo "  filter        Filter BAM reads based on kinetic signals and other metrics."
  echo "  map           Map filtered BAM to a reference genome."
  echo "  detect_dz     Detect dZ-DNA reads using MLP and SVM classifiers."
  echo "  classify      Perform single-nucleotide Z/A classification on BAM."
  echo "  tax_source    Identify the taxonomic source of dZ-DNA reads.
  --input_type  Mandatory input type (mapped/unmapped). If 'unmapped', the BAM will be mapped to a reference genome before zfreq analysis.
  --ref_species  Optional reference species file for read-ref analysis"
  echo "  convert       Convert BAM with ZP tags into a tab-delimited table."
  echo "  zfreq         Calculate Z/(A+Z) at each reference A/T base."
  echo
  echo "Options (for all subcommands):"
  echo "  -i, --input       Input file (BAM or FASTA)"
  echo "  -o, --output      Output file"
  echo "  -r, --ref         Reference FASTA file (required for mapping)"
  echo "  -m, --model       Model for MLP or SVM operations"
  echo "  -t, --threads     Number of threads to use"
  echo "  -h, --help        Display this help message"
}

# Exit if no arguments are provided
if [ $# -eq 0 ]; then
  usage
  exit 1
fi

# Parse the subcommand
subcommand=$1
shift

# Initialize variables for options
input=""
output=""
ref=""
model=""
threads=1

# Parse options
while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      input="$2"
      shift 2
      ;;
    -o|--output)
      output="$2"
      shift 2
      ;;
    -r|--ref)
      ref="$2"
      shift 2
      ;;
    -m|--model)
      model="$2"
      shift 2
      ;;
    -t|--threads)
      threads="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --read_ref)
      read_ref="$2"
      shift 2
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

# Check if input file is provided
if [ -z "$input" ]; then
  echo "Error: Input file is required."
  usage
  exit 1
fi

# Perform subcommands based on user input
case $subcommand in
  filter)
    if [ -z "$output" ]; then
      echo "Error: Output file is required for filtering."
      usage
      exit 1
    fi
    echo "Filtering BAM reads..."
    # Apply filters to BAM reads
    python py/filter_bam.py -b "$input" -o "$output"
    ;;

  zfreq)
    if [ "$input_type" = "unmapped" ]; then
      echo "Mapping filtered BAM to reference..."
      # Map filtered BAM to reference
      pbmm2 align "$ref" "$input" "$output" --preset CCS --sort
      input="$output"
    fi
    if [ -z "$input_type" ] || [ -z "$output" ]; then
      echo "Error: Input type and output files are required for zfreq."
      usage
      exit 1
    fi
    
    if [ "$input_type" = "unmapped" ] && [ -z "$ref" ]; then
      echo "Error: Reference file is required for unmapped input type."
      usage
      exit 1
    fi
      usage
      exit 1
    fi
    echo "Mapping filtered BAM to reference..."
    # Map filtered BAM to reference
    pbmm2 align "$ref" "$input" "$output" --preset CCS --sort
    ;;

  detect_dz)
    if [ -z "$model" ]; then
      echo "Error: Model is required for dZ-DNA detection."
      usage
      exit 1
    fi
    echo "Detecting dZ-DNA reads..."
    # Run MLP module to detect dZ-DNA reads
    z-calling-base/z-calling-base -k 21 -t "$threads" "$input" "$model" "$output"
    # Run SVM classifier on reads containing ZP tag
    z-calling-read/z-calling-read predict "${output%bam}${model}.bam" "${output%bam}${model}.SVM.tsv" model/svm/k21ReadClassifier 500 100
    ;;

  classify)
    if [ -z "$model" ] || [ -z "$output" ]; then
      echo "Error: Model and output files are required for classification."
      usage
      exit 1
    fi
    echo "Classifying Z/A bases..."
    # Run MLP module for single-nucleotide Z/A classification
    z-calling-base/z-calling-base -k 11 -t "$threads" "$input" "$model" "$output"
    ;;

  tax_source)
    ref_species=""
read_ref=""
    if [ -z "$ref" ] || [ -z "$output" ]; then
      echo "Error: Reference and output files are required for taxonomic source analysis."
      usage
      exit 1
    fi
    echo "Identifying taxonomic source of dZ-DNA reads..."
    # Convert BAM to FASTA
    samtools fasta -@ "$threads" "$input" > "${input%bam}fasta"
    # Identify taxonomic source using minimap2
    minimap2 -t "$threads" -x map-hifi "$ref" "${input%bam}fasta" | awk '!a[$1]++' | awk '($4-$3/$2)>0.4 {match($0, /dv:f:([0-9.]+)/, a); if (a[1]<0.05) {print $1"\t"$6}}' > read_RefContig.tsv
    python py/SVM_LR_analysis.py --SVM "${input%bam}${model}.SVM.tsv" --output "$output" --read_ref read_RefContig.tsv ${ref_species:+--ref_species $ref_species}
    ;;

  convert)
    if [ -z "$output" ]; then
      echo "Error: Output file is required for BAM to table conversion."
      usage
      exit 1
    fi
    echo "Converting BAM to tab-delimited table..."
    # Convert BAM to tab-delimited table with ZP tags
    z-bam2txt/z-bam2txt -zp -t "$threads" "$input" - > "$output"
    ;;

  zfreq)
    if [ -z "$output" ]; then
      echo "Error: Output file is required for Z frequency calculation."
      usage
      exit 1
    fi
    echo "Calculating Z frequency..."
    # Calculate Z/(A+Z) for each reference A/T base
    z-freq/zfreq "$input" "$output"
    ;;

  *)
    echo "Error: Unknown subcommand: $subcommand"
    usage
    exit 1
    ;;
esac

# End of script
exit 0
