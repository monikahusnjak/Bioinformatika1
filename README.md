# Bioinformatics Variant Calling Pipeline

This project implements a simple variant calling pipeline in C++ using Minimap2, Samtools mpileup, a custom mpileup parser, and a custom variant caller.

The pipeline was tested on two datasets:

* Lambda genome
* E-coli genome


## Requirements

The following software must be installed:

sudo apt update
sudo apt install g++ samtools minimap2

## Clone Repository

git clone https://github.com/monikahusnjak/Bioinformatics-1.git
cd Bioinformatics-1


## Large Files

The file:

ecoli/data/ecoli_simulated_reads.fasta

is not included in the repository because of its size.

Download it from:

Releases -> Bioinformatics Project Dataset

and place it in:

ecoli/data/ecoli_simulated_reads.fasta

before running the E. coli pipeline.


## Running the E. coli Pipeline

cd ecoli
chmod +x run.sh
./run.sh

The script performs:

1. Read alignment using Minimap2
2. SAM to BAM conversion
3. BAM sorting
4. BAM indexing
5. mpileup generation
6. mpileup parsing
7. Variant calling

Output files are generated automatically.


## Running the Lambda Pipeline

cd lambda
chmod +x run.sh
./run.sh

The same workflow is applied to the lambda dataset.


## Pipeline Overview

Reference genome
        ↓
Read alignment (Minimap2)
        ↓
       SAM
        ↓
       BAM
        ↓
    sorted BAM
        ↓
     mpileup
        ↓
  mpileup parser
        ↓
candidate variants
        ↓
 variant caller
        ↓  
 final variants
