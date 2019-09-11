##GWAS data

I downloaded some example data from dbGaP. These data were for running through methods only and do not represent a valid GWAS comparison. Here are the commands I used do download these data.

```sh
# download ped and map for 1008 TB cases
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83397/suppl/GSE83397_project1008s.ped.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83397/suppl/GSE83397_project1008s.map.gz
# download 1000G SNP chip data to use as controls
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz.tbi
#download some additional data from Chinese population run on the same chip as TB cases to use as controls
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69664/suppl/GSE69664_GPL19864_processed_data.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69664/suppl/GSE69664_GPL20166_processed_data.txt.gz
```

I then converted all data to plink binary format, mapped to reference strand, and merged keeping only SNPs in the intersection of all data sets.
