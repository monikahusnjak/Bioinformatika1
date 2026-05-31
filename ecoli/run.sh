#!/bin/bash

set -e

echo "Running E. coli pipeline..."

echo "1. Aligning reads with minimap2..."
minimap2 -x map-ont -a --eqx data/ecoli.fasta data/ecoli_simulated_reads.fasta > data/aln.sam

echo "2. Converting SAM to BAM..."
samtools view -Sb data/aln.sam > data/aln.bam

echo "3. Sorting BAM..."
samtools sort data/aln.bam -o data/aln.sorted.bam

echo "4. Indexing BAM..."
samtools index data/aln.sorted.bam

echo "5. Generating mpileup..."
samtools mpileup -f data/ecoli.fasta data/aln.sorted.bam > pileup.txt

echo "6. Compiling mpileup parser..."
g++ -std=c++17 stmpileup_parser.cpp -o stmpileup_parser

echo "7. Running mpileup parser..."
./stmpileup_parser

echo "8. Compiling variant caller..."
g++ -std=c++17 ecoli_variants.cpp -o ecoli_variants

echo "9. Running variant caller..."
./ecoli_variants

echo "E. coli pipeline finished."
