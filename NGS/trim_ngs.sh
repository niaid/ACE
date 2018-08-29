#!/bin/bash

## adapter - ktrim=r means trim on the right
## ref=/path/to/bbduk/adapters.fa, overwrite=t saying it's ok to overwrite any files named the same as the out files, threads=3 use 3 processors
bbduk.sh in=SRR2057563_1.fastq in2=SRR2057563_2.fastq out=SRR2057563_trimmed.1.fastq out2=SRR2057563_trimmed.2.fastq ktrim=r ref=/opt/anaconda/opt/bbmap-38.22-0/resources/adapters.fa overwrite=t threads=3 



## quality - can do at the same time as adapter, just did this separately for clarity
## qtrim=rl - trim quality on right and left
## trimq=15 - trim threshold score 15
bbduk.sh in=SRR2057563_trimmed.1.fastq in2=SRR2057563_trimmed.2.fastq out=SRR2057563_trimmed.quality.1.fastq out2=SRR2057563_trimmed.quality.2.fastq qtrim=rl trimq=15 overwrite=t threads=3 





## primer to trim amplicon data - we didn't do this in the class
bbduk.sh in=22057_S2_R1_subsample.fastq in2=22057_S2_R2_subsample.fastq out=22057_S2_R1_subsample_trimmed.fastq out2=22057_S2_R1_subsample_trimmed.fastq ktrim=l k=75 mink=18 ref=/home/ace/ace_workshop/ngs_basics_extra/primers.fa copyundefined=t minlen=60 overwrite=t threads=3
