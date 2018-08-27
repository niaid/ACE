I. RNA-seq introduction

A. RNA-seq goals

  1. Functional studies

     *a) Differential gene expression*

  2. SNP analysis

  3. Transcriptome assembly

  4. Novel gene finding

  5. Splice variant analysis

B. Background - sequencing

  1. Illumina

  2. Library prep

  3. Paired vs single end

C. Sequencing files

  1. Fastq

II. RNA-seq workflow 

A. Introduction 

  1. Type of sequencing

  2. Experimental design

     _a) block design_

     _b) coverage considerations_

  3. Eukaryote vs prokaryote

     *a) Tools for eukaryote*
       
       - Hisat2 - string tie - ballgown

       - STAR - DESeq2

     _b) Prokaryote_

       - Bowtie2 - HTseq - DESeq2

B. basic steps

  1. QC

     _a) QC trimming? trimmomatic_

     _b) Adaptor trimming_

     _c) use python script for example_

  2. Mapping

     *a)  Alignment output*
        
        - SAM and BAM

        - Show reads mapped in IGV

        - Introduction to samtools

        - Flagstats

  3. Counting - use python script

C. DESeq2 introduction

  1. Distribution assumptions

  2. Count based

  3. Statistics used

D. DESeq2 workshop

  1. Importing experiment design

  2. Importing read matrix

  3. Exploratory

      _a) PCA_

      _b) Count matrix heatmaps_

  4. Differential expression

      _a) Extract significant up and down_

  5. Example showing how to separate out all the plasmids

  6. Examples on how to add annotation

  7. Briefly introduce pathway analysis, GO terms.  Mariam will go over
    this later
