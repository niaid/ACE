#!/bin/bash

# Author: Amitava Roy
# Purpose: The document is a script which utilizes the while-loop structure of bash
#          and introduces the floating point math in bash scripting. 
#          Note that bash does not support floating point arithmetic.
#          However, we can pipe other executable, like python, to carry out our math.
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
# bash while_loop.sh

COUNTER=0
while [  $COUNTER -lt 10 ]; 
do
    echo print $COUNTER/3. | python
    let COUNTER=COUNTER+1 
done

