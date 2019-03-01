#!/bin/bash

### This script is written to be used on Locus 
# To run it elsewhere, you will need to install the following:
# 1. GNU parallel https://www.gnu.org/software/parallel/
# 2. docopts for argument parsing (found in the class repo)
# 2. readlength.sh from BBMap suite from JGI https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/
#
# As written, it's expected that this script, docopts, and readlength.sh are in the same directory.
# To run elsewhere, you should change the paths to wherever the tools are installed on your computer.
# 

## directory where this script, docopts, and readlength.sh are located
dir=$(dirname "$0")

## path to where docopts is located - change if needed to wherever you have docopts
docopts=${dir}/docopts



eval "$(${docopts} -b -h - : "$@" <<EOF
Usage: $0 [options] <listfile>
Options:
   --r1 <r1>             suffix for R1 or unpaired; listfile has R1 or unpaired filenames [default: _R1_001.fastq]
   --r2 <r2>             suffix for R2; default is same as r1 with R2 instead
   -u, --unpaired        reads are unpaired; listfile has all filenames
   -o <outfile>     	 output file [default: seqinfo.txt]
   -k, --keep       	 keep intermediate files
   -t <nprocs>      	 number of processors [default: 4]
   -s, --separatepairs   produce table with R1 R2 separated
EOF
)"



if [ -z ${r2} ]; then
    r2=${r1/R1/R2}
    echo ${r2}
fi

## Load the parallel module; comment out if running 
module load parallel

if [[ ${unpaired} == 'true'  ]]; then
    cat ${listfile} | parallel -j ${nprocs} ${dir}/readlength.sh in={} out={="s/.+\/|${r1}//g"=}.seqinfo.txt
elif [[ ${separatepairs} == 'true'  ]]; then
    cat ${listfile} | parallel -j ${nprocs} ${dir}/readlength.sh in={} out={="s/.+\/|${r1}//g"=}.seqinfo.txt
    cat ${listfile} | parallel -j ${nprocs} ${dir}/readlength.sh in={="s/${r1}//"=}${r2} out={="s/${r1}//"=}${r2}.seqinfo.txt
else
    cat ${listfile} | parallel -j ${nprocs} ${dir}/readlength.sh in={} in2={="s/${r1}//"=}${r2} out={="s/.+\/|${r1}//g"=}.seqinfo.txt
fi


if [[ ${separatepairs} == 'true'  ]]; then
    echo -e "Sample\tReads\tBases.R1\tMax.R1\tMin.R1\tAvg.R1\tMedian.R1\tMode.R1\tStd_Dev.R1\tBases.R2\tMax.R2\tMin.R2\tAvg.R2\tMedian.R2\tMode.R2\tStd_Dev.R2" >${outfile}
    cat ${listfile} | parallel -q perl -ane 'if ($_ =~ /Reads/) { $s = $ARGV; $s =~ s/\.seqinfo\.txt//; print "$s"; } last if $. > 8; print "\t$F[1]"; END { print "\n"; }' {="s/.+\/|${r1}//g"=}.seqinfo.txt  >${outfile}.1
    cat ${listfile} | parallel -q perl -ane 'if ($_ =~ /Reads/) { $s = $ARGV; $s =~ s/\.seqinfo\.txt//; print "$s"; } last if $. > 8; print "\t$F[1]"; END { print "\n"; }' {="s/${r1}/${r2}/g"=}.seqinfo.txt | cut -f 3-  >${outfile}.2
    paste ${outfile}.1 ${outfile}.2 >>${outfile}
else
    echo -e "Sample\tReads\tBases\tMax\tMin\tAvg\tMedian\tMode\tStd_Dev" >${outfile}

    cat ${listfile} | parallel -q perl -ane 'if ($_ =~ /Reads/) { $s = $ARGV; $s =~ s/\.seqinfo\.txt//; print "$s"; } last if $. > 8; print "\t$F[1]"; END { print "\n"; }' {="s/.+\/|${r1}//g"=}.seqinfo.txt  >>${outfile}
fi

if [ ${keep} == 'false' ]; then
    cat ${listfile} | parallel rm {="s/.+\/|${r1}//g"=}.seqinfo.txt
    if [ ${separatepairs} == 'true' ]; then
	cat ${listfile} | parallel rm {="s/${r1}/${r2}/g"=}.seqinfo.txt
	rm ${outfile}.2 ${outfile}.1
    fi
fi


