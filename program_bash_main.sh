#!/bin/bash

#check if necessary arguments are introduced
if [ $# -ne 2 ]; then
 echo "The program needs to arguments, genome file and reads file (in this order)";
 exit;
fi
GENOME=$1
FILE=$2
export GENOME
export FILE

# Make temp directory for files produced
if [ ! -d assign4_temp ]; then
  mkdir assign4_temp;
fi

FIRSTCHAR=$(cat $FILE | head -c 1);
# If first character is > run fasta script
if [ $FIRSTCHAR = ">" ]; then
    printf "\nInput sequence identified as fasta type\n"
    source program_bash_fasta.sh;
fi

# If first character is @ run fastq script
if [ $FIRSTCHAR = "@" ]; then
    printf "\nInput sequence identified as fastq type\n"
    source program_bash_fastq.sh;
fi



