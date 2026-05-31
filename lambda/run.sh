#!/bin/bash

set -e

echo "Running lambda pipeline..."

echo "1. Aligning reads with minimap2..."
minimap2 -x map-ont -a --eqx data/lambda.fasta data/lambda_simulated_reads.fasta > data/aln.sam

echo "2. Converting SAM to BAM..."
samtools view -Sb data/aln.sam > data/aln.bam

echo "3. Sorting BAM..."
samtools sort data/aln.bam -o data/aln.sorted.bam

echo "4. Indexing BAM..."
samtools index data/aln.sorted.bam

echo "5. Generating mpileup..."
samtools mpileup -f data/lambda.fasta data/aln.sorted.bam > pileup.txt

echo "6. Compiling mpileup parser..."
g++ -std=c++17 stmpileup_parser.cpp -o stmpileup_parser

echo "7. Running mpileup parser..."
./stmpileup_parser

echo "8. Compiling variant caller..."
g++ -std=c++17 lambda_variants.cpp -o lambda_variants

echo "9. Running variant caller..."
./lambda_variants

echo "Lambda pipeline finished."
