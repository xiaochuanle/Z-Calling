# Z-Calling

## Brief Introduction

Z-Calling is a machine-learning based tool for (1) distinguishing dZ-DNA molecules from ordinary DNAs and (2) calling Z bases in A/Z-coexisting DNAs. Z-Calling also has implemented a pipeline for detecting taxonomic or sequence sources of dZ-DNAs in mixed datasets. In multiple tested datasets, Z-Calling has faithfully identified souces of dZ-DNAs without false positive discory. And its Z base calling module has achieved AUCs ranging from 0.9422 to 0.9550 and F1 scores ranging from 0.85-0.92 across all tested datasets.

## Building from Scratch

### Preparing the Python Environment

Create a virtual environment using Conda. The Python scripts require numpy (version 20.0 or higher) and pytorch (version 2.0 or higher) with CUDA 11.8 support.

```bash
conda create -n Z-Calling python=3.11
conda activate Z-Calling
pip install numpy torch==2.0.1
```

### Building the C++ Program

Z-Calling was tested and runed in **NVIDIA GeForce RTX 3090**,  ensure you have a **GPU** and **CUDA Toolkit 11.8** installed.  Download **libtorch 2.0.1** if it's not already included in your Python environment. This C++ program is compiled using g++-11.2 on Ubuntu 22.04.

**If you are not familiar about how to install CUDA Toolkit 11.8, here is a example for set up CUDA Toolkit 11.8 in ubuntu 22.04 x86_64 system**

```bash
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-ubuntu2204.pin
sudo mv cuda-ubuntu2204.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/11.8.0/local_installers/cuda-repo-ubuntu2204-11-8-local_11.8.0-520.61.05-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu2204-11-8-local_11.8.0-520.61.05-1_amd64.deb
sudo cp /var/cuda-repo-ubuntu2204-11-8-local/cuda-*-keyring.gpg /usr/share/keyrings/
sudo apt-get update
sudo apt-get -y install cuda
```

**Install the following packages before building the program:**

1. boost

2. zlib

**And the these projects are already included in `3rdparty/`**

4. [argparse](https://github.com/p-ranav/argparse "argparse"): Argument Parser for Modern C++
5. [htslib](https://github.com/samtools/htslib "htslib"): An implementation of a unified C library for accessing common file formats
6. [libsvm](https://github.com/cjlin1/libsvm "libsvm"): An efficient software for SVM classification and regression.


```bash
cd z-calling-base
mkdir build && cd build
cmake -DCMAKE_PREFIX_PATH=`python -c 'import torch;print(torch.utils.cmake_prefix_path)'` .. # Determine the cmake path # if you haven`t set up the python environment, you should directy include libtorch path here.
make -j

cd z-calling-read
mkdir build && cd build
cmake ..
make -j

cd z-bam2txt
mkdir build && cd build
cmake ..
make -j

cd z-freq
mkdir build && cd build
cmake ..
make -j

cd z-seq
mkdir build && cd build
cmake ..
make -j
```

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
cd z-calling-base
./z-calling-base -k Kmer_size -t Threads_num Bam_path Model_path Output_bam_path

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
cd z-calling-read
./z-calling-read predict Bam_path Result_path Model_path length_threshold section_num

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
cd z-bam2txt
./z-bam2txt -zp -t Threads_num Bam_path - > Result_path

Positional arguments:
  Threads_num           Number of CPU threads (default: 16)
  Bam_path              Path to the Z-base called BAM file
  Result_path           Path to result tsv file
```

### Covert Z-Called BAM into six-alphabet fasta file
```bash
cd z-seq
./z-seq Bam_path Result_path

Positional arguments:
  Bam_path              Path to the Z-base called BAM file
  Result_path           Path to result fasta file
```

### Calculate Z/(A+Z) at each reference A/T bases for a mapped bam
```bash
cd z-freq
./zfreq input_tsv output_tsv

Positional arguments:
  input_tsv           Path to the converted tab delimited table file
  output_tsv          Path to the output table file
```
