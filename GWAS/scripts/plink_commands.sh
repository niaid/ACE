#!/bin/bash

##start by checking counts of individuals and SNPs using "wc -l " and R for Ca/Co counts

##initial data filtering
base_dir="/etc/ace-data/ace-training/karlins.eric/GWAS_lecture_Sept2019"
#in a clean directory
cp -r $base_dir .
cd GWAS_lecture_Sept2019
ls -r *
#it's mostly empty directories

##start by checking counts of individuals and SNPs using "wc -l " and R for Ca/Co counts

plink --bfile plink_start/TB_GWAS --mind 0.05 --geno 0.05 --maf 0.05 --hwe 0.0001 --allow-no-sex --make-bed --out qc1/TB_GWAS

#check counts of snps and individuals after filtering

#QC filtered, right?  Let's try associatin test

plink --bfile qc1/TB_GWAS --assoc fisher --adjust --allow-no-sex --out assoc1/TB_GWAS

##switch to R script and make manhattan and qq plots

##run pca
plink --bfile qc1/TB_GWAS --pca 6 header --allow-no-sex --out pca1/TB_GWAS

## switch to R script and plot pca

## subset ancestry matched samples
plink --bfile qc1/TB_GWAS --keep ancestry_filt1/samples_to_keep.txt --make-bed --out ancestry_filt1/TB_GWAS

## ok , we ancestry matched let's run the association test again

plink --bfile ancestry_filt1/TB_GWAS --assoc fisher --adjust --allow-no-sex --out assoc2/TB_GWAS

#can make manhattan and qqplots again or just look at lambda value in log file (stdout from plink)

# run pca again

plink --bfile ancestry_filt1/TB_GWAS --pca 6 header --allow-no-sex --out pca2/TB_GWAS

## go back to R and plot pca

# subset from the optmatch

plink --bfile ancestry_filt1/TB_GWAS --keep ancestry_filt2/samples_to_keep.txt --make-bed --out ancestry_filt2/TB_GWAS

# run assoc test

plink --bfile ancestry_filt2/TB_GWAS --assoc fisher --adjust --allow-no-sex --out assoc3/TB_GWAS

#go to R and make final GWAS plots

