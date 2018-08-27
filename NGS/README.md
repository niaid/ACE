# NGS Basics

#### A. Short History
- Sanger
- NGS - Illumina (SbS), Ion Torrent (semiconductor), 454 (pyrosequencing) - defunct 
- Third Generation - PacBio, Oxford nanopore

#### B. File formats
- common formats and converting from/to different formats
- FASTQ/FASTA/bam/bed
- demultiplexing

#### C. QC
- FASTQC/multiQC - gauge quality
  - common problems - primer/adapter dimer; R2 poor quality
  - https://blog.horizondiscovery.com/diagnostics/the-5-ngs-qc-metrics-you-should-know
- trimming
  - adapters/primers and quality
  - when should you do it?
  - common tools, example using bbduk

#### D. Mapping 
- Why do you map reads to a reference? 
  - reference guided assembly (DNA-seq or RNA-seq - transcript assembly), variant calling, removing host background/contaminants, gene expression (RNA-seq)
- Reference sequence - where do you get it, assess its quality, etc
- short read mapping - bwa, bowtie2
  - common parameters - identity, coverage, gaps
  - example using bowtie2
- long read mapping - minimap2, MashMap
- RNA-seq - HiSAT, STAR
- algorithms? - should we discuss this? if this is supposed to be more hands on and we only have 1.5 hrs, not sure if we will have time.  Can include information in materials
- examining the alignment
  - coverage, paired-end "happiness", % reads matching
  - bedtools, samtools, picardtools
  - IGV, Tablet

#### E. Practical

- to be done all the way through the talk?
- we can start with some small sample fastq files and go through all the steps up to mapping against a small reference (bacteria, virus?)













