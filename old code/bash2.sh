#!/bin/bash
#merge sam
# directory assign4_output contains all the split .fasta sequences
# creates merged sam file from all the separate sam files in assign4_output
for file in assign4_output/*.sam; do grep -v '^@' $file >> merged-file.sam ;done
