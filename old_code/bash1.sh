#!/bin/bash
# directory assign4_output contains all the split .fasta sequences
# sacch.fasta - complete saccharomyces genome
bwa index -p sacch -a bwtsw sacch.fasta;
for f in assign4_output/*.fasta; do bwa aln -t 4 sacch $f > $f.bwa| bwa samse sacch $f.bwa $f > $f.sam ;done;
