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
7.Compute how many reads have been aligned (using Linux tools)
"""
from itertools import product
import re


class Program():
    def __init__(self, input_file_name):
        # Open input file, create list with lines from this file
        with open(input_file_name, 'rt') as input_file:
            input_file_list = input_file.readlines()
            input_file_list = [line.strip() for line in input_file_list]
            file_extension = input_file_name.split('.')[1]

        # Transform Fastq input into lists of list containing heads and reads respectively
        if input_file_list[0].startswith('@') and input_file_list[2] == '+':
            head = [item.replace("@", ">") for item in input_file_list[::4]]
            read = [item for item in input_file_list[1::4]]
            qual = [item for item in input_file_list[3::4]]

        # Transform Fasta input into lists of list containing heads and reads respectively
        elif input_file_list[0].startswith('>') and bool(re.match('^[ATCGN]+$', input_file_list[1].upper())):
            # Create list of lists with different repeated lines from FASTA file
            head = [item for item in input_file_list[::2]]
            read = [item for item in input_file_list[1::2]]

        # Save variables to use in the rest of the program
        self.head = head
        self.read = read



    # Counts all the 3-mers occuring in the fasta/fastq file
    def mer_abundance(self):
        # Create list of all the possible 3-mers from ATCGN
        possible_mers = [''.join(c) for c in product('ATCGN', repeat=3)]
        # Intiate dictionary with possible mers as keys and values as 0 (for counting of mers)
        mer_dict = dict(zip(possible_mers, [0] * len(possible_mers)))
        # Split each sequence line every 3 characters and add to 3-mer counts in the dictionary
        for line in self.read:
            for mer in re.findall('...', line):
                mer_dict[mer] += 1
        # Return dictionary containing all 3-mer counts
        return mer_dict


    # Splits each header+sequence into a seperate fasta file in the output directory

    # COMMENTED IT OUT FOR NOW SO IT DOESN'T PRODUCE BUNCH OF FILES IF YOU TRY TO RUN THE PROGRAM.
    '''
    def split_to_files(self):
        for i, header in enumerate(self.head):
            # Each file will be named fasta_seq and contain corresponding index to the header or seq list

            #with open(f"output/fasta_seq{i}", 'wt') as output_file:
            #    output_file.write(f'{header}\n{self.read[i]}\n')
                
        # Print number of files produced (== number of reads)
        print(f"\nNumber of fasta files produced: {len(self.head)}")
    '''

    # Hard trims all the sequences 20 nucleotides from the right,
    # If the sequence is shorter or equal than 20 nucleotides, append 0 instead for later processing
    def hard_trim(self):
        new_read = []
        for seq in self.read:
            if len(seq) > 20:
                new_seq = seq[:-20]
            else:
                new_seq = 0
            new_read.append(new_seq)

        # Save list of trimmed reads as a class variable.
        self.new_read = new_read


if __name__ == "__main__":
    # Load sample fastq file
    input_file_name = "sample.fastq"
    # Initiate class
    prog = Program(input_file_name)
    # Perform mer_abundance function from class
    prog.mer_abundance()
    # Perform splitting sequence from fasta/fastq file into separate transformed fasta files
    #prog.split_to_files()
    # Hard trim all the sequences 20 nucleotides from the right
    prog.hard_trim()
