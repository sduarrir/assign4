#!/bin/bash
# directory assign4_output contains all the split .fasta sequences
# sacch.fasta - complete saccharomyces genome
bwa index -p sacch -a bwtsw sacch.fasta;
for file in assign4_output/*.fasta; do bwa aln -t 4 sacch $file > $file.bwa| bwa samse sacch $file.bwa $file > $file.sam ;done