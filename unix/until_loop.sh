#!/bin/bash

# Author: Amitava Roy
# Purpose: The document is a script which utilizes the until-loop structure of bash
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
# bash until_loop.sh

COUNTER=20
while [  $COUNTER -ge 10 ]; 
do
    echo COUNTER $COUNTER
    let COUNTER-=1 
done

