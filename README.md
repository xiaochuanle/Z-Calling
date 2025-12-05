# Z-Calling

## Brief Introduction

Z-Calling is a machine-learning toolkit designed to detect and call 2,6-diaminopurine (Z-base) modifications using PacBio HiFi reads.

## Core Features

- **Molecule Classification:** Distinguishes dZ-DNA molecules from canonical (A) DNA.
- **Base Calling:** Precisely identifies Z bases within hybrid A/Z sequences.
- **Source Detection:** A pipeline to identify taxonomic sources of dZ-DNAs in mixed datasets with zero false positive discovery in benchmarking.

## Performance

Across multiple datasets, the Z base calling module achieves **AUCs of 0.94–0.96** and **F1 scores of 0.85–0.92**.

## Computational Requirements & Installation

### System Requirements
Z-Calling is designed to run on standard Linux environments and does **not** require GPU acceleration.

* **Operating System:** Linux
* **Memory (RAM):** ≥ 32 GB recommended (Peak memory usage is approx. 31 GB during base calling)
* **Processor:** Standard Multi-core CPU (Benchmarked on AMD EPYC 7402 24-Core Processor)
* **Disk Space:** Sufficient space for input BAMs and output files

### Building from Scratch
Z-Calling is designed to run on CPU. The installation process uses Conda (or Mamba) to handle all dependencies automatically, including HTSlib, PyTorch, and the necessary C++ compilers. This eliminates the need for manual library compilation.

#### 1. Download this repository.
```base
git clone https://github.com/xiaochuanle/Z-Calling.git
cd Z-Calling
```

#### 2. Set up the Environment
We strongly recommend using Mamba for faster environment resolution, though standard Conda will also work. The environment.yml file in the root directory of the project.
Run the following commands to create the environment:

```bash
# Option A: Using Mamba (Recommended - Faster)
mamba env create -f environment.yml

# Option B: Using Standard Conda (Slower)
conda env create -f environment.yml
```
Once the environment is created, activate it and test if torch is successfully installed:
```bash
conda activate Z-Calling
python -c "import torch; print(torch.__version__)"  # Successful installation will print '2.2.2+cpu'
```

#### 3. Building the Program
```bash
mkdir build
cd build
cmake ..
make -j
```
Once compilation is complete, the executables (such as z-calling-base and z-bam2txt) will be located in the build/ directory.
Optionally: Add Z-Calling/build to your PATH.

## Z-Calling Usage

### Filters were applied to filter BAM reads
#### 1. kinetic signals (fi, ri, fp, rp tags) are required. Otherwise Z-Calling will quit.
#### 2. Reads with (Median_PW / median_IPD < 0.3) or (Median_PW / median_IPD > 99percentile) were filtered.
#### 3. The default minimum pass number for forward (fn) and reverse (rn) passs is 3 (fn>=3 and rn>=3).

```bash
python py/filter_bam.py -b RawBAM -o FilteredBAM

Positional arguments:
  RawBAM            Input bam could be either mapped or unmapped
  FilteredBAM       Output path to the filtered BAM file    

Optional: Mapping filtered BAM to a reference fasta file
pbmm2 align REF FilteredBAM FilteredMappingBAM --preset CCS --sort
```

### Z-calling for base and read
#### Running MLP module
```bash
build/z-calling-base -k Kmer_size -t Threads_num Bam_path Model_path Output_bam_path

Positional arguments:
  Kmer_size              K-mer size for Z calling (21 when using k21-full-ZA model, and 11 when using k11-mixed-AZ model)
  Threads_num            Number of CPU threads (default: 16)
  Bam_path               Path to the input BAM file
  Model_path             Path to the trained model
                         Running MLP module with k21-full-ZA model. dZ-DNA read detection: the model was trained on full-dA/dZ datasets and intended to be used for dZ-DNA read classification only.
                         Running MLP module with k11-mixed-AZ model. Single-nucleotide Z/A classification: a ZP tag will be added to reads that records the Z probability for each A/T base in the sequence.
  Output_bam_path        Path to the result BAM file
```
#### Running SVM classifier on reads containing ZP tag.
```bash
build/z-calling-read predict Bam_path Result_path Model_path length_threshold section_num

Positional arguments:
  Bam_path              Path to the Z-base called BAM file
  Result_path           Path to result tsv file
  Model_path            Path to the trained model 
  length_threshold      The minimum number of bases (default: 500)
  section_num           The number of intervals (default: 100)
```

#### Taxonomic source of dZ-DNA reads
```bash
samtools fasta -@8 FilteredBAM > Filteredfasta
minimap2 -t 8 -x map-hifi $REF Filteredfasta | awk '!a[$1]++' | awk '($4-$3/$2)>0.4 {match($0, /dv:f:([0-9.]+)/, a); if (a[1]<0.05) {print $1"\t"$6}}' > read_RefContig.tsv
python py/SVM_LR_analysis.py --SVM z-calling-read-Result_path --read_ref read_RefContig.tsv [ --ref_species RefContig_Spcecies.tsv ] --output SVM.LR.tsv
```

#### Optional: convert ZP tag bam into tab delimited table file. Per A/T Z probabilities were reported in six columns (unmapped) or eight column format (mapped, corresponding reference coordinate also reported).
```bash
build/z-bam2txt -zp -t Threads_num Bam_path - > Result_path

Positional arguments:
  Threads_num           Number of CPU threads (default: 16)
  Bam_path              Path to the Z-base called BAM file
  Result_path           Path to result tsv file
```

### Covert Z-Called BAM into six-alphabet fasta file
```bash
build/z-seq Bam_path Result_path

Positional arguments:
  Bam_path              Path to the Z-base called BAM file
  Result_path           Path to result fasta file
```

### Calculate Z/(A+Z) at each reference A/T bases for a mapped bam
```bash
build/zfreq input_tsv output_tsv

Positional arguments:
  input_tsv           Path to the converted tab delimited table file
  output_tsv          Path to the output table file
```
