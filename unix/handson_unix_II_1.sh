#!/bin/bash

# Author: Amitava Roy    
# Purpose: The document contains some of the basic commands of UNIX
#          examples of their usage.
# Date -
# Creation : August 12, 2018
# Last modification : August 12, 2018
#
# Usage -
# The document is for reference purpose only. If needed it can be used as
# bash handson_unix_II_1.sh  
# Some of the materials has been used from the book 
# Unix and Perl Primer for Biologists (Version 3.1.2 — October 2016) - Keith Bradnam & Ian Korf

# In this excercise we will need two terminals open. In one terminal we will go through this
# file. Inthe second one we will create some other files.

# Let us create our first script. Open a file hello.sh and type the following lines
##!/bin/bash ! This line tells the OS that it is a bash shell script. Other shell languages are csh, tcsh, ksh ...
# 
# echo "Hello World"

# Now try running
bash hello.sh

# Try running as an executable
./hello.sh

# Programs in Unix need permission to be run. We will normally always have to type the following for any script that
# we create:
chmod u+x hello.sh

# Try running again as an executable. ./ means the current directory.
./hello.sh

# This would use the chmod to add executable permissions (+x) to the file called ‘hello.sh’ (the ‘u’ means add this
# permission to just you, the user). Without it, your script won’t run. The above chmod command would add executable 
# permissions (+x) to the file called ‘hello.sh’ (the ‘u’ means add this permission to just you, the user).
# The chmod command can also modify read and write permissions for files, and change any of the three sets of
# permissions (read, write, execute) at the level of ‘user’, ‘group’, and ‘other’. You probably won’t need to know any
# more about the chmod command other than you need to use it to make scripts executable.
# Let us look at the chmod commands
man chmod

# Now let us try to execute the hello.sh again.
hello.sh

# It will not run as the hello.sh is in a directory where the UNIX is not instructed to look for an executable.
# Let us open the .profile file in your home directory and add the current directory to the PATH.
# If you do not have the .profile file, let us creat onthe UNIX is not instructed to look for an executable.
# Let us open the .profile file in your home directory and add the current directory to the PATH.
# If you do not have the .profile file, let us creat one and add
# PATH=$PATH":your current directory"

source ~/.profile
hello.sh

# Let us write a script, organize.sh to organize files in the directory. The goal of he script will be to
# 1) Create two directories "Sequences" and "Scripts"
# 2) Copy all files with extension .fasta to the directory "Sequences"
# 3) Copy all the files with extension .sh to the directopry "Scripts"

# Now that we have an example of if-then-else-fi in organize.sh, we will examine the different loop structures in bash 
# shell. Let us examine the file for_loop.sh.

# Next we will examine more traditional for loop strucutre along with the arithmetica capability in bash with 
# the script for_loop_c.sh.

# Not that bash does not support floating point arithmetic. If we need to carry out a floating poitn operation
# we can use other executables withing a bash script, like python, to carry out the operation. An example
# is given in while_loop.sh. The script also examines the while-do-done structure of bash.

# Finally we look at the until-do-done strucutre in bashin the script until_loop.sh, 
# which is a while-do-done like but opposite in function.
