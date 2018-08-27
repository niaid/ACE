#!/bin/bash

# Author: Amitava Roy    
# Purpose: The document contains some of the basic commands of UNIX
#          examples of their usage.
# Date -
# Creation : August 11, 2018
# Last modification : August 11, 2018
#
# Usage -
# The document is for reference purpose only. If needed it can be used as
# bash handson_unix_I_2.sh  
# Some of the materials has been used from the book 
# Unix and Perl Primer for Biologists (Version 3.1.2 — October 2016) - Keith Bradnam & Ian Korf

# You might already know that FASTA files (used frequently in bioinformatics) have a simple format: one
# header line which must start with a ‘>’ character, followed by a DNA or protein sequence on subsequent lines. To
# find only those header lines in a FASTA file, we can use grep, which just requires you specify a pattern to search for,
# and one or more files to search:

grep ">" b_burgdorferi_cp9.fasta

# This will produce lots of output which will flood past your screen. If you ever want to stop a program running in Unix,
# you can type Control+C (this sends an interrupt signal which should stop most Unix programs). The grep command
# has many different command-line options (type man grep to see them all), and one common option is to get grep
# to show lines that don’t match your input pattern. You can do this with the -v option and in this example we are
# seeing just the sequence part of the FASTA file.

grep -v ">" b_burgdorferi_cp9.fasta

# By now, you might be getting a bit fed up of waiting for the grep command to finish, or you might want a cleaner
# way of controlling things without having to reach for Ctrl-C. Ideally, you might want to look at the output from any
# command in a controlled manner, i.e. you might want to use a Unix program like less to view the output.
# This is very easy to do in Unix, you can send the output from any command to any other Unix program (as long as
# the second program accepts input of some sort). We do this by using what is known as a pipe. This is implemented
# using the ‘|’ character (which is a character which always seems to be on different keys depending on the keyboard
# that you are using). Think of the pipe as simply connecting two Unix programs. In this next example we send the
# output from grep down a pipe to the less program. Let’s imagine that we just want to see lines in the input file
# which contain the pattern “ATGTGA” (stop codon):

grep "TGA" b_burgdorferi_cp9.fasta | less

# Notice that you still have control of your output as you are now in the less program. If you press the forward slash
# (/) key in less , you can then specify a search pattern. Type TGA after the slash and press enter. The less
# program will highlight the location of these matches on each line. Note that grep matches patterns on a per line
# basis. So if one line ended TG and the next line started A, then grep would not find it.

# Sometimes we do not want to use less to see all of the output from a command like grep. We might just want to
# see a few lines to get a feeling for what the output looks like, or just check that our program (or Unix command) is
# working properly. There are two useful Unix commands for doing this: head and tail. These commands show (by
# default) the first or last 10 lines of a file (though it is easy to specify more or fewer lines of output). So now, let’s look
# for another pattern which might be in all the sequence files in the directory. If we didn’t know whether the
# DNA/protein sequence in a FASTA files was in upper-case or lower-case letters, then we could use the -i option of
# grep which ‘ignores’ case when searching. Note - for simple searches you do not have to use inverted commas:

grep -i TGA b_burgdorferi_cp9.fasta | head
grep -i TGA b_burgdorferi_cp9.fasta | tail

# A concept that is supported by many Unix programs and also by most programming languages (including Python) is
# that of using regular expressions. These allow you to specify search patterns which are quite complex and really
# help restrict the huge amount of data that you might be searching for to some very specific lines of output. E.g. you
# might want to find lines that start with an ‘AAA’ and lines finish with ‘TTT’:

grep "^AAA" b_burgdorferi_cp9.fasta
grep "TTT$" b_burgdorferi_cp9.fasta

# The ^ character is a special character that tells grep to only match a pattern if it occurs at the start of a line. 
# Similarly, the $ tells grep to match patterns that occur at the end of the line. Now let us try to find lines that  
# start with an ‘AAA’ and finish with ‘TTT’

grep "^AAA.*TTT$" b_burgdorferi_cp9.fasta

# Now let us try to find lines that and finish with ‘TTT’, but which have at least two GTT in the middle 

grep "^AAA.*GTTGTT.*TTT$" b_burgdorferi_cp9.fasta

# The . and * characters are also special characters that form part of the regular expression. Try to understand
# how the following patterns all differ. Try using each of these these patterns with grep against any one of the
# sequence files. Can you predict which of the five patterns will generate the most matches?
# ACGT
# AC.GT
# AC*GT
# AC.*GT

# The asterisk in a regular expression is similar to, but NOT the same, as the other asterisks that we
# have seen so far. An asterisk in a regular expression means: ‘match zero or more of the preceding
# character or pattern’.
# AC.GT searches for AC[A/T/G/C]GT
# AC*GT searches for AGT, ACGT, ACCGT, ACCCGT etc.
# AC.*GT searches for AC[any number of any character]GT

# Rather than showing you the lines that match a certain pattern, grep can also just give you a count of how many
# lines match. This is one of the frequently used grep options. Running grep -c simply counts how many lines
# match the specified pattern. It doesn’t show you the lines themselves, just a number:

grep -c "AC.*GT" b_burgdorferi_cp9.fasta

# Unix has a very powerful command called sed that is capable of performing a variety of text manipulations. 
# Let’s assume that you want to change the way the FASTA header looks: 

head -n 1 b_burgdorferi_cp9.fasta

# You want to change the word "cp9" to "copy 9" and make all the lower case "b" to upper case "B".
# We can request multiple substitutions with each preceded by the option "-e"
# Option "g" at the end of each substitutions indicate that all instances of the query will be substituted.
# Otherwise only the first instance is substituted.
# Note that this doesn’t actually change the contents of the file, it just changes the screen output from the previous
# command in the pipe.

head -n 1 b_burgdorferi_cp9.fasta | sed -e 's/cp9/copy 9/' -e 's/b/B/g'

# We can send the output to a file

sed -e 's/cp9/copy 9/' -e 's/b/B/g' b_burgdorferi_cp9.fasta > temp.fasta

# This step introduces a new concept. Up till now we have sent the output of any command to the screen (this is the
# default behavior of Unix commands), or through a pipe to another program. Sometimes you just want to redirect the
# output into an actual file, and that is what the ">" symbol is doing, it acts as one of three redirection operators in
# Unix.
