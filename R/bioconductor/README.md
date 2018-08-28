# Bioconductor Exercises

Most of these exercises are taken from past workshops given at the BioC conferences which are put on the [Bioconductor team](http://bioconductor.org/) annually.

Each topic has 2 (4) files (except for ggplot2).

1. *topic.nb.html* - R notebook which has the information for importing/creating the dataset and the exercises.  
   - *topic.Rmd* is the markdown file used to create the html - it can also be downloaded from the html page in your browser using the dropdown menu on the top right (so, you can also just download the html)
2. *topic_solution.nb.html* - R notebook that also contains *the solutions* to the exercises (we don't have this for ggplot2, but I can add soon)

### R notebooks

- [**Summarized Experiment**](summarized_experiment.Rmd): exercises that use the SummarizedExperiment data structure.  The first set is introductory; the second set combines with some alignment data and is a little more advanced.
- [**Introduction to ggplot2**](ggplot2-introduction.Rmd): (Very) short introduction to ggplot2.  Uses [data.csv](data.csv) as input dataset.
- [**BAM files**](bamfiles.Rmd): use GenomicAlignments bioconductor package to read in and examine a bam file.
- [**Genomic Ranges**](genomiciranges.Rmd): constructing and manipulating a GRanges data structure.

### Data files

- We use the data from 2 bioconductor packages:

  - [pasillaBamSubset](http://bioconductor.org/packages/release/data/experiment/html/pasillaBamSubset.html)
  - [TxDb.Dmelanogaster.UCSC.dm3.ensGene](http://bioconductor.org/packages/release/data/annotation/html/TxDb.Dmelanogaster.UCSC.dm3.ensGene.html)

- For the workshop, we put gzipped files in the `ace_bioconductor` folder, but you can also install directly from bioconductor (see links above)

- To install from the gzipped file within R:

  ```R
  install.packages("/path/to/packagefile.gz", repos=NULL)
  ## for the ACE Uganda workshops the path and commands should be:
  install.packages("/home/ace/ace_workshop/ace_bioconductor/pasillaBamSubset_0.18.0.tar.gz", repos=NULL)
  install.packages("/home/ace/ace_workshop/ace_bioconductor/TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2.tar.gz", repos=NULL)
  ```
