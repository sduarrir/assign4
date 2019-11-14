# -*- coding: utf-8 -*-
"""
@author: Jakub Widawski
@email: jdwidawski@gmail.com
@python: Python3.7
            ---     TASK    ---
Produce two scripts (Bash/Python) that, given an input FASTA/FASTQ file,
implement the following pipeline:
1.Counts the abundance of each possible 3-mers (histogram)
2.Splits the input in FASTA/FASTQ into files of only 1 sequence each
3.Hard-trims (from the right) all sequences from all files 20nt
4.Using the genome mapping tool BWA and the reference genome of the ScaromiceCerevisiae(any strain will do), aligns each of the files producing its corresponding SAM file.
5.Merge all SAM files ignoring headers (using Linux tools)
6.Sorts the SAM file by chromosome and position
7.Compute how many read have been aligned (using Linux tools)
"""
from itertools import product
import re
import collections
from Bio import SeqIO
import os

class Program():
    def __init__(self, input_file_name):
        # Open input file, create list with lines from this file
        with open(input_file_name, 'rt') as input_file:
            input_file_list = input_file.readlines()[:8000] # only use first 8000 lines for testing so its fast
            input_file_list = [line.strip() for line in input_file_list]

        # Transform Fastq input into lists of list containing head and read respectively
        if input_file_list[0].startswith('@'):
            filetype = "fastq"
            heads = [item for item in input_file_list[::4]]
            reads = [item.upper() for item in input_file_list[1::4]]
            quals = [item for item in input_file_list[3::4]]
            seq_dict = {}
            for header,seq,q in zip(heads, reads, quals):
                seq_dict[header] = [seq,q]


        # Transform Fasta input into lists of list containing head and read respectively
        elif input_file_list[0].startswith('>'):
            filetype = "fasta"
            records = SeqIO.parse(input_file_name, "fasta")
            seq_dict = {}
            for record in records:
                header = record.description
                seq = ''.join(record.seq.upper())
                q = "0"*len(record.seq)
                seq_dict[header] = [seq,q]

        else: #we should have into account non fasta or fastaq (printing an error message)
            raise ValueError("Program just works with fasta and fastq format")
            #dont know if value error is the most accurate type of error

        # Save variables to use in the rest of the program

        self.heads = list(seq_dict.keys())
        self.reads = [x[0] for x in seq_dict.values()]
        self.seq_dict = seq_dict
        self.filetype = filetype


    # Counts all the 3-mers occuring in the fasta/fastq file
    def mer_abundance(self):
        # Create list of all the possible 3-mers from ATCGN
        possible_mers = [''.join(c) for c in product('ATCG', repeat=3)]
        # Intiate dictionary with possible mers as keys and values as 0 (for counting of mers)
        mer_dict = dict(zip(possible_mers, [0] * len(possible_mers)))
        # Split each sequence line every 3 characters and add to 3-mer counts in the dictionary
        for line in self.reads:
            for mer in re.findall('...', line):
                try:
                    mer_dict[mer] += 1
                except KeyError:
                    pass

        sorted_x = sorted(mer_dict.items(), key=lambda kv: kv[1])
        mer_dict = collections.OrderedDict(sorted_x)
        print("\n 3-mer counts in the input sequence:\n")
        [print("\t", key, ' ==> ', value) for key, value in mer_dict.items()]
        # Return dictionary containing all 3-mer counts
        #return mer_dict


    # Splits each header+sequence into a seperate fasta file in the assign4_temp directory
    # COMMENTED IT OUT FOR NOW SO IT DOESN'T PRODUCE BUNCH OF FILES IF YOU TRY TO RUN THE PROGRAM.
    def split_to_files(self):
        if self.filetype == "fasta":
            for i, header in enumerate(self.heads):
                # Each file will be named fasta_seq and contain corresponding index to the header or seq list
                with open(f"out_py/fasta_seq{i}.fasta", 'wt') as output_file:
                    new_seq = self.seq_dict[header][0]
                    output_file.write(f"{header}\n{new_seq}\n")

            # Print number of files produced (== number of reads)
            print(f"\nNumber of fasta files produced: {len(self.heads)}")

        elif self.filetype == "fastq":
            for i, header in enumerate(self.heads):
                # Each file will be named fasta_seq and contain corresponding index to the header or seq list
                with open(f"out_py/fastq_seq{i}.fastq", 'wt') as output_file:
                    new_seq = self.seq_dict[header][0]
                    qual = self.seq_dict[header][1]
                    output_file.write(f"{header}\n{new_seq}\n+\n{qual}")


    '''
    # Hard trims all the sequences 20 nucleotides from the right,
    # If the sequence is shorter or equal than 20 nucleotides, append 0 instead for later processing
    def hard_trim(self):
        file_list = [f for f in os.listdir('out_py') if os.path.isfile(os.path.join('out_py', f))]
        for filename in file_list:
            if filename.endswith(".fasta"):
                print(filename)
                if len(seq) > 20:
                    new_seq = seq[:-20]
                else:
                    new_seq = 0
                new_reads.append(new_seq)
            elif filename.endswith(".fastq"):
                with open(filename, 'r+') as input_file:
                    input_file_list = input_file.readlines()
                    print(input_file_list)
                    seq = input_file_list[1]
                    if len(seq) > 20:
                        new_seq = seq[:-20]
                    print(seq)
                    break
    '''




if __name__ == "__main__":
    # Load sample fastq file
    input_file_name = "sequence.fasta"
    # Initiate class
    prog = Program(input_file_name)
    # Perform mer_abundance function from class
    mer_dictionary = prog.mer_abundance()

    # Perform splitting sequence from fasta/fastq file into separate transformed fasta files
    prog.split_to_files()
    # Hard trim all the sequences 20 nucleotides from the right
    #prog.hard_trim()