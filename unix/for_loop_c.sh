#!/bin/bash

# Author: Amitava Roy
# Purpose: The document is a script which utilizes the for-loop structure of bash.
#          The script creates a C-like for loop and introduces the math in bash scripting
#
# Some notes about arithimatic in bash
# let expression - Make a variable equal to an expression.
# expr expression - print out the result of the expression.
# $(( expression )) - Return the result of the expression.
# ${#var} - Return the length of the variable var.
#
# Date -
# Creation : August 12, 2018
# Last modification : August 12, 2018
#
# Usage -
# bash for_loop_c.sh

for i in `seq 1 10`; 
do
   let j=i+1
   k=$(( $i * $i ))
   l=${#k}
   echo $i $j $k $l
done

