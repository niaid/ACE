#!/bin/bash

# author: 		mariam quinones
# purpose: 		This script will allow students to run a set of minimal steps required to perform variant analysis 
#				of a trio (mom, dad and daughter) dataset. 
# workflow: 	1) map reads to the human genome; 2) Process bam files; 3) Call variants using gatk; and 4) Annotate variants using snpEff
# dependencies:	GATK, snpEff, samtools, IGV, 


#############################################################################
#			 Step 0: Create the working directory 							#
#############################################################################

# original data source: http://figshare.com/articles/6_files_with_1GB_per_file/106340 
# Corpas Family Trio Exome Data. Manuel Corpas. figshare.
# http://dx.doi.org/10.6084/m9.figshare.106340
# Retrieved 14:13, Aug 06, 2013 (GMT)

# Note: The original data files were modified (reduced to 500k reads each) 
#	 	to be used during the workshop without the need for intense computation
# 		The modified files are available for download here: 
#		https://s3.amazonaws.com/ace-uganda/exome.corpas.tar.gz

mkdir ~/Documents/variants_run
cp -R ~/ace_workshop/share_files/share_corpas/* ~/Documents/variants_run/

# create directory for input reads
cd ~/Documents/variants_run
mkdir map_out
mkdir input_reads/
cd input_reads/

mv ../*clean5.33.fq .


#############################################################################
#			 Step 1: Map the reads to the human genome chromosome 16	    #
#############################################################################

#Note:		We will only use chr 16 of the entire human genome (version hg19) 
#			for this exercise but you should always map to the entire genome

# As good practice, run fastqc and multiqc to examine quality.
fastqc *
multiqc .

# Create aliases for the genome file and the output directory
humanref=~/Documents/variants_run/human_g1k_v37.chr16.fasta
destdir=~/Documents/variants_run/map_out

# Run this loop to map reads using bwa and then create a bam file using samtools
# 	Important note: We have already created an index for the reference genome, 
# 	otherwise you would first need to run the following lines to generate all dictionary and indeces:
#	bwa index ~/Documents/variants_run/human_g1k_v37.chr16.fasta
# 	samtools faidx human_g1k_v37.chr16.fasta
#	picard CreateSequenceDictionary R=human_g1k_v37.chr16.fasta O=human_g1k_v37.chr16.dict

for fname in *_1.clean5.33.fq
do
	base=${fname%_*1*}
	bwa mem $humanref -R "@RG\\tID:${base}\\tSM:${base}\\tPL:ILLUMINA" -v 2 -t 4 "${base}_1.clean5.33.fq" "${base}_2.clean5.33.fq" > "$destdir/${base}-pe.5.33.sam"
	samtools view -b -T $humanref -o $destdir/${base}-pe.5.33.bam $destdir/${base}-pe.5.33.sam
	samtools sort $destdir/${base}-pe.5.33.bam > $destdir/${base}-pe.sorted.5.33.bam
	samtools index $destdir/${base}-pe.sorted.5.33.bam
	samtools idxstats $destdir/${base}-pe.sorted.5.33.bam
	echo "$base"
done




#############################################################################
#		Step 2: Preprocess the bam file to recalibrate base quality         #
#############################################################################

# give to gatk, multiple known variant files (vcf) as appropriate to train the algorithm
# for this demo we will only use 1000 genomes but it is recommended to use dbsnps or others.


#dbsnp=~/Documents/variants_run/dbsnp_138.b37.vcf
gold=~/Documents/variants_run/Mills_and_1000G_gold_standard.indels.b37.vcf

# optional: limit the variant search to certain locations by providing a 'bed' file
regions=~/Documents/variants_run/S03723314_Regions_chr16.fix.bed

cd $destdir
ls *-pe.sorted.5.33.bam | sed 's/-.*//g' | sort | uniq > ID.txt

file="$destdir/ID.txt"
  while read line
  do
	picard MarkDuplicates INPUT=$destdir/$line-pe.sorted.5.33.bam OUTPUT=$destdir/$line.chr16.dedup.bam AS=true CREATE_INDEX=true M=$destdir/$line.chr16.metrics.txt
	gatk BaseRecalibrator -I $destdir/$line.chr16.dedup.bam -R $humanref --known-sites $gold --output $destdir/$line.recal_data.table -L 16
	gatk ApplyBQSR -I $destdir/$line.chr16.dedup.bam -R $humanref -L 16 --bqsr-recal-file $destdir/$line.recal_data.table --output $destdir/$line.chr16.dedup.recal.bam
	echo "$line"
done <"$file"

# verify bam files using browser IGV

#############################################################################
#				Step 3: Variant calling using GATK						    #
#############################################################################

# Use GATK Haplotype caller to derive likelihoods of genotypes and creates gVCF
file="$destdir/ID.txt"
  while read line
  do
	gatk HaplotypeCaller -R $humanref -ERC GVCF --input $destdir/$line.chr16.dedup.recal.bam --output $destdir/$line.chr16.g.vcf
	echo "$line"
done <"$file"

ls *chr16.g.vcf > outputlist.txt

# combine g.vcfs and reate a database
gatk GenomicsDBImport \
    -V mom.chr16.g.vcf \
    -V daughter.chr16.g.vcf \
    -V dad.chr16.g.vcf \
    --genomicsdb-workspace-path corpas_database \
    --intervals 16
    
# use GenotypeGVCFs to create a merged vcf file and apply filters if desired ; note the gendb:// prefix to the database input directory path.
gatk GenotypeGVCFs \
    -R $humanref \
    -V gendb://corpas_database \
    -O merged_corpas.chr16.g.vcf 
    

#############################################################################
#				Step 4: Variant annotation using SnpEFF					    #
#############################################################################

# Note: We have already downloaded the hg19 database using the line of code right below
#java -Xmx4G -jar /hpcdata/bcbb/quinones/software/snpEff_Aug2018/snpEff/snpEff.jar download hg19

# Run SnpEff with one line which generates the annotated vcv and a report in html format
java -Xmx4G -jar /hpcdata/bcbb/quinones/software/snpEff_Aug2018/snpEff/snpEff.jar hg19 merged_corpas.chr16.g.vcf > merged_corpas.chr16.g.annotated.vcf






##########################################################################################################################################################
#	Additional potential steps to refine variant calling.   These are demonstrated as if using GATK 3; they will need to be slightly modifed to run with gatk 4
##########################################################################################################################################################

#!!!!!!!!!!!! DO NOT RUN!!!!!!!!!!!*************
# 3) GATK VQSR Variant recalibration- WILL NOT DO
		##if you have a cohort, we will not run this for this exercise####################
# 		gatk VariantRecalibrator \
# 		  -R $humanref \
# 		  -input merged_corpas.chr16.g.vcf \
# 		  -L 16\
# 		  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ../hapmap_3.3.b37.vcf  \
# 		  -resource:omni,known=false,training=true,truth=false,prior=12.0 ../1000G_omni2.5.b37.vcf \
# 		 -resource:1000G,known=false,training=true,truth=false,prior=10.0 ../1000G_phase1.snps.high_confidence.b37.vcf \
# 		 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ../dbsnp_138.b37.vcf  \
# 		 -an QD -an FS -an SOR  \
# 		 -mode SNP \
# 		 -recalFile SNP.recal \
# 		 -tranchesFile SNP.tranches \
# 		 -rscriptFile SNP.plots.R
# 
# 		gatk ApplyRecalibration \
# 		 -R ../human_g1k_v37.fasta  \
# 		 -input  trio_corpas_gatk.vcf \
# 		 -mode SNP \
# 		 -recalFile SNP.recal \
# 		 -tranchesFile SNP.tranches \
# 		 -o trio_corpas_gatk.recal.SNP.vcf \
# 
# 
# 		do the same for indel
# 		gatk VariantRecalibrator \
# 		 -R ../human_g1k_v37.fasta \
# 		 -input trio_corpas_gatk.recal.SNP.vcf \
# 		 -L ../S03723314_Regions_chr16.fix.bed \
# 		 -resource:mills,known=false,training=true,truth=true,prior=12.0 ../Mills_and_1000G_gold_standard.indels.b37.vcf \
# 		 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ../dbsnp_138.b37.vcf \
# 		 -an QD -an FS -an SOR -an FS -an ReadPosRankSum \
# 		 -mode INDEL \
# 		 -recalFile INDEL.recal \
# 		 -tranchesFile INDEL.tranches \
# 		 -rscriptFile INDEL.plots.R
# 
# 		gatk ApplyRecalibration \
# 		 -R ../human_g1k_v37.fasta  \
# 		 -input  trio_corpas_gatk.recal.SNP.vcf \
# 		 -mode INDEL \
# 		 -recalFile INDEL.recal \
# 		 -tranchesFile INDEL.tranches \
# 		 -o trio_corpas_gatk.recal.SNP_INDEL_RECAL.vcf \
# 		 --ts_filter_level 99.0
# 
# 		############################################################################################################
# 		END
# 		##########################################
# 		**** If VQSR is not applicable ( <30 samples)****########
# 
# 		extract SNPs and apply hard filter, do the same for indels
# 		extracting SNPs
# 		gatk SelectVariants \
# 			-R ../human_g1k_v37.fasta  \
# 			-L ../S03723314_Regions_chr16.fix.bed \
# 			-V sample.vcf \
# 			-selectType SNP \
# 			-o sample_SNPs.vcf  
# 		applying hard filter
# 		
#		gatk VariantFiltration \
# 			-R ../human_g1k_v37.fasta  \
# 			-L ../S03723314_Regions_chr16.fix.bed \
# 			-V sample_SNPs.vcf \
# 			--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
# 			--filterName "my_snp_filter" \
# 			-o sample.filtered_SNPs.vcf 
# 
# 		extracting INDELs
# 		gatk SelectVariants \
# 			-R ///human_g1k_v37.fasta  \
# 			-V sample.vcf \
# 			-selectType INDEL \
# 			-o sample.indels.vcf 
# 
# 		applying hard filters to indels
# 		gatk VariantFiltration \
# 			-R ../human_g1k_v37.fasta \
# 			-V sample.indels.vcf \
# 			--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
# 			--filterName "my_indel_filter" \
# 			-o sample.filtered_indels.vcf 
# 
# 		merge two files
# 
# 		gatk CombineVariants \
# 		 -L ../S03723314_Regions_chr16.fix.bed \
# 		 --variant:indel sample.filtered_indels.vcf \
# 		 --variant:snps sample.filtered_SNPs.vcf \
# 		 -R ../human_g1k_v37.fasta  \
# 		 -o sample.filtered_SNP_indels.vcf \
# 		 -genotypeMergeOptions PRIORITIZE \
# 		 -priority snps,indel

##################################################




    