## Commands used in NGS Basics Presentation

We have 2 datasets in our directory:

- Shotgun data
  - SRR2057563_1.fastq and SRR2057563_2.fastq
- 16S amplicon data - fastq files beginning with 220 and ending in subsample.fastq
  - 220*_subsample.fastq

**File formats**

To convert from fastq to fasta, we use `fastq_to_fasta` from the FASTX toolkit

```bash
fastq_to_fasta -n -v -i 22057_S2_R1_subsample.fastq -o 22057_S2_R1_subsample.fasta
```

**FastQC**

We run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](http://multiqc.info/) to get some nice graphs of our data quality.

We will make some directories to store our output.

```bash
mkdir fastqc1
mkdir fastqc2
```

Run fastqc and multiqc on the SRR shotgun data

```
fastqc SRR*.fastq -o ./fastqc1
cd ./fastqc1
multiqc .
```

Run fastqc on the 220 amplicon data

```bash
fastqc 220*.fastq -o ./fastqc2
cd ./fastqc2
multiqc .
```

The output will be called *multiqc_report.html* and will be in folders *fastqc1* and *fastqc2*

**Trimming**

- [trim_ngs.sh](trim_ngs.sh) - script with trimming commands

- rerun fastqc to compare the trimmed reads to the original reads

```bash
fastqc SRR2057563_trimmed.quality*.fastq -o ./fastqc1
cd fastqc1
multiqc -f .
```

How does the trimmed quality differ from the untrimmed?

**Mapping**

We will map our trimmed reads against a viral genome.  The genome is this one: [Cowpox virus](https://www.ncbi.nlm.nih.gov/assembly/GCF_000839185.1/) with refseq accession GCF_000839185.1

We use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) for the mapping.

1. build the cowpox bowtie2 index

   We give the FASTA filename of the cowpox genome and the name of the index, which we will name "cowpox."  This will create files with extension *.bt2

```bash
bowtie2-build ./GCF_000839185.1_ViralProj14174_genomic.fna cowpox
```

2. inspect bowtie2 cowpox index

```bash
bowtie2-inspect -s cowpox
```

3. Map against the cowpox index

```bash
## simple mapping
bowtie2 -x cowpox -1 SRR2057563_trimmed.quality.1.fastq -2 SRR2057563_trimmed.quality.2.fastq -S SRR2057563_cowpox.bowtie2.sam

## more sensitive mapping
## see the bowtie2 manual for what the parameters mean! http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
bowtie2 -D 20 -R 3 -p 3 -N 1 -L 25 -i S,1,0.50 -t -x cowpox -1 SRR2057563_1.fastq -2 SRR2057563_2.fastq -S SRR2057563_cowpox.bowtie2.sam --no-unal 2>&1 | tee SRR2057563_cowpox.bowtie2.log
```

This produces a sam file *SRR2057563_cowpox.bowtie2.sam*

4. sort the sam file

bowtie2 produces an unsorted sam file - the coordinates of the alignments in the reference is jumbled.  A lot of programs that use sam files (including igv) require a sorted sam file

```bash
samtools sort SRR2057563_cowpox.bowtie2.sam -o SRR2057563_cowpox.bowtie2.sorted.sam
```



5. Look at alignments in igv

- import viral genome in **Genomes -> Load Genome from File...** Navigate to the directory with the genome.  On the ACE laptop, the file is */home/ace/ace_workshop/ace_ngs_basics_data/GCF_000839185.1_ViralProj14174_genomic.fna*

- import gff (annotations) with **File -> Load from File...** Choose the gff file */home/ace/ace_workshop/ace_ngs_basics_data/GCF_000839185.1_ViralProj14174_genomic.gff*

- import the alignment with **File -> Load from File..**. Choose *SRR2057563_cowpox.bowtie2.sorted.sam*

- right click left side of alignment track

  - color by insert size pair orientation
6. Picard tools
  - get insert size metric for igv
- ```bash
  picard CollectInsertSizeMetrics I=SRR2057563_cowpox.bowtie2.sorted.sam O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5
  ```

- In igv - right click left side of alignment track -> set insert size options...

- Group alignments by pair orientation

  - <https://broadinstitute.github.io/picard/explain-flags.html>

- Want to filter out unmapped pairs

  - samtools to filter

```bash
samtools view -h -F 8 SRR2057563_cowpox.bowtie2.sorted.sam -o SRR2057563_cowpox.properlypaired.sam
```

7. Remove human reads from original files - for submission to NCBI, for example

```bash
bowtie2 -k 1 -D 20 -R 3 -p 3 -N 1 -L 25 -i S,1,0.50 -t -x Homo_sapiens.GRCh38.dna.chromosome.1 -1 SRR2057563_1.fastq -2 SRR2057563_2.fastq -S SRR2057563_hg38.bowtie2.sam --no-unal --un-conc SRR2057563_unmapped.fastq 2>&1 | tee SRR2057563_hg38.bowtie2.log
```

