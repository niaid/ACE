#!/usr/bin/env bash

mkdir fastqc_output

## run fastqc
fastqc SRR2057563_*.fastq -o ./fastqc_output
cd ./fastqc_output
multiqc .

## inspect bowtie2 human index
bowtie2-inspect -s GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index

## map against human
bowtie2 -k 1 -D 20 -R 3 -p 3 -N 1 -L 25 -i S,1,0.50 -t -x GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 SRR2057563_1.fastq -2 SRR2057563_2.fastq -S SRR2057563_hg38.bowtie2.sam --no-unal --un-conc SRR2057563_unmapped.fastq 2>&1 | tee SRR2057563_hg38.bowtie2.log

## build cowpox bowtie2 index
bowtie2-build ./GCF_000839185.1_ViralProj14174_genomic.fna cowpox

## map against bos mutus


## map against cowpox
bowtie2 -D 20 -R 3 -p 3 -N 1 -L 25 -i S,1,0.50 -t -x cowpox -1 SRR2057563_unmapped.1.fastq -2 SRR2057563_unmapped.2.fastq -S SRR2057563_cowpox.bowtie2.sam --no-unal 2>&1 | tee SRR2057563_cowpox.bowtie2.log


## samtools sort
samtools sort SRR2057563_cowpox.bowtie2.sam -o SRR2057563_cowpox.bowtie2.sorted.sam
