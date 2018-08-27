#!/bin/bash

# Author: Amitava Roy
# Purpose: The document is a script which utilizes the for-loop structure of bash.
#          The script lists the subdirectories in the current directory.
#
# Date -
# Creation : August 12, 2018
# Last modification : August 12, 2018
#
# Usage -
# bash for_loop.sh

for i in $( ls ); 
do
   if [ -d "$i" ]; then
     echo subdirectory: $i
   fi
done

