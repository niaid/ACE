#!/usr/bin/env bash

## Class datasets can be copied to your own directory on Locus:
cp -r /hpcdata/bcbb/poorani/NGSclass/qc_data .
cd ./qc_data



#### QC #####
module load fastqc/0.11.5-Java-1.8.0_45
module load multiqc/1.5-Python-3.6.0

#### FastQC
mkdir fastqc_output

## fastqc help
fastqc -h

## run fastqc
fastqc *.fastq -o ./fastqc_output
cd ./fastqc_output

## run multiqc to generate summary report
## useful parameter to add is --cl-config "read_count_multiplier: 1"
## this will give the total number of reads in each file instead of giving the count in millions
multiqc .



#### Trimming
module load bbmap/38.41

## Move back up to the qc_data directory
cd ..

## BBDuk help https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
bbduk.sh -h

## primer trimming on the left
bbduk.sh in=poor_16S_S1_R1_001.fastq in2=poor_16S_S1_R2_001.fastq out=poor_16S_S1_R1_trimmed.fastq out2=poor_16S_S1_R2_trimmed.fastq ktrim=l k=18 mink=4 ref=./primers.fa copyundefined=t overwrite=t threads=10 tbo


## adapter trimming on the right 
bbduk.sh in=SRR2057563_1.fastq in2=SRR2057563_2.fastq out=SRR2057563_trimmed.1.fastq out2=SRR2057563_trimmed.2.fastq ktrim=r k=27 mink=4 hdist=1 ref=../bbmap/resources/adapters.fa minlen=10 overwrite=t threads=10


## quality - can do at the same time as adapter, just did this separately for clarity
bbduk.sh in=SRR2057563_trimmed.1.fastq in2=SRR2057563_trimmed.2.fastq out=SRR2057563_trimmed.quality.1.fastq out2=SRR2057563_trimmed.quality.2.fastq qtrim=rl trimq=15 maq=2 minlen=60 overwrite=t threads=10



## SRA data is from from https://www.ebi.ac.uk/ena/data/view/SRR2057563
## data was sequenced from sample collected from cow lesions containing cowpox virus; so may have some cow DNA and the cowpox virus


#### Mapping ######
module load bowtie2/2.3.4.1

## See the bowtie2 help http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
bowtie2 --help


## map against bos taurus to remove background host reads
## we used --very-fast in class, but probably better to omit this parameter or use --very-sensitive
## there is trade-off between speed and sensitivity
bowtie2 -k 1 -p 10 --very-fast -t -x /hpcdata/bcbb/poorani/NGSclass/bowtie2/bostaurus -1 SRR2057563_trimmed.quality.1.fastq -2 SRR2057563_trimmed.quality.2.fastq -S SRR2057563_bostaurus.bowtie2.sam --no-unal --un-conc SRR2057563_unmapped.fastq

## bos taurus index was built from GCF_002263795.1_ARS-UCD1.2_genomic.fna
## can use my index (in command above)
## OR can build your own index; this make take some time
## To build your own, modify the cowpox index command below to instead build from the Bos taurus genome

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
