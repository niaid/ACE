#!/usr/bin/env bash


### add bbduk to your path
export PATH=$PATH:/classhome/bbmap

#### QC #####
module load fastqc/0.11.5-Java-1.8.0_45
module load multiqc/1.5-Python-3.6.0

#### FastQC
mkdir fastqc_output

## run fastqc
fastqc SRR2057563_*.fastq -o ./fastqc_output
cd ./fastqc_output
multiqc .

#### Trimming


## adapter
bbduk.sh in=SRR2057563_1.fastq in2=SRR2057563_2.fastq out=SRR2057563_trimmed.1.fastq out2=SRR2057563_trimmed.2.fastq ktrim=r k=27 mink=4 hdist=1 ref=/classhome/bbmap/resources/adapters.fa minlen=10 overwrite=t threads=10


## primer
bbduk.sh in=22057_S2_R1_subsample.fastq in2=22057_S2_R2_subsample.fastq out=22057_S2_R1_subsample_trimmed.fastq out2=22057_S2_R1_subsample_trimmed.fastq ktrim=l k=75 mink=18 ref=/home/ace/ace_workshop/primers.fa copyundefined=t minlen=60 overwrite=t threads=10

## quality - can do at the same time as adapter, just did this separately for clarity
bbduk.sh in=SRR2057563_trimmed.1.fastq in2=SRR2057563_trimmed.2.fastq out=SRR2057563_trimmed.quality.1.fastq out2=SRR2057563_trimmed.quality.2.fastq qtrim=rl trimq=15 maq=2 minlen=60 overwrite=t threads=10




#### Mapping ######
module load bowtie2/2.3.4.1

## map against bos mutus
bowtie2 -k 1 -p 10 --very-fast -t -x /classhome/bowtie2/bostaurus -1 SRR2057563_1.fastq -2 SRR2057563_2.fastq -S SRR2057563_bostaurus.bowtie2.sam --no-unal --un-conc SRR2057563_unmapped.fastq

## build cowpox bowtie2 index
bowtie2-build ./GCF_000839185.1_ViralProj14174_genomic.fna cowpox

## inspect bowtie2 index
bowtie2-inspect -s cowpox



## map against cowpox
bowtie2 -t -x cowpox -1 SRR2057563_unmapped.1.fastq -2 SRR2057563_unmapped.2.fastq -S SRR2057563_cowpox.bowtie2.sam --no-unal 2>&1 | tee SRR2057563_cowpox.bowtie2.log


#### Samtools
module load samtools/1.9-goolf-1.7.20

## samtools sort
samtools sort SRR2057563_cowpox.bowtie2.sam -o SRR2057563_cowpox.bowtie2.sorted.sam
