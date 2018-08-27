#!/bin/bash

# Author: Amitava Roy
# Purpose: The document is a script which will create directories Seqeunces and Scripts if they do not exist.
#          An examples of if-then-else-fi inbuilt structure in bash script.
# Date -
# Creation : August 12, 2018
# Last modification : August 12, 2018
#
# Usage -
# bash organize.sh

if [ -d "Sequences" ]; then
    echo "Directory alread exists. Not creating a directory."
else
    mkdir Sequences
fi
cp *.fasta Sequences/

if [ -d "Scripts" ]; then
    echo "Directory alread exists. Not creating a directory."
else
    mkdir Scripts
fi
cp *.sh Scripts/

