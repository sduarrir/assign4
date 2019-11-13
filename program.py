"""
            ---     TASK    ---
Produce two scripts (Bash/Python) that, given an input FASTA/FASTQ file, implement the following pipeline:
"""
import sys
import os
from Bio import SeqIO #would allow us to process sequences
##FUNCTIONS


# Function to open input file and check if correct format
def open_input(filename): #we provide the file
    try:

    except FileNotFoundError:
        print(filename, "not found. check for any typos!")
        exit(1)
    line = inputfile.readline()
    if line[0] == "@":
        fileformat = "fastq"
    elif line[0] == ">":
        fileformat = "fasta"
    else:
        print(filename, "format is not fastaq or fasta.") #if format is not fasta or fastaq it will print an error and exit
        exit(1)
    inputfile.seek(0)
    return (inputfile, fileformat) #returns the file and its format (fasta o fastaq)


#open the files (and check if format is correct
if len(sys.argv) !=3:
    print("The program need the following arguments:\n\t--input [INPUTFILE]\n\t--ouput [FILE]  If file exists, will be overwrited!!!\n\t--operation [OPERATION] (what do you want the program to print\nSome operations will need conditions.\n\trc (reverse complement) does not need conditions\n\ttrim needs: --trim-right [POSITION] and --trim-left [POSITION]\n\tadaptor-removal needs: --adaptor [SEQUENCE]")
    exit(1)
genome = sys.argv[1]
reads = sys.argv[2]

##1.Counts the abundance of each possible 3-mers (histogram)



##2.Splits the input in FASTA/FASTQ into files of only 1 sequence each
##3.Hard-trims (from the right) all sequences from all files 20nt
##4.Using the genome mapping tool BWA and the reference genome of the ScaromiceCerevisiae(any strain will do), aligns each of the files producing its corresponding SAM file.

##5.Merge all SAM files ignoring headers (using Linux tools)
os.system("samtools merge merge.sam *")

##6.Sorts the SAM file by chromosome and position
os.system("samtools sort merge.sam sorted_merged.sam")
os.system("samtools index sorted_merged.sam ind_sorted_merged.sam")

##7.Compute how many reads have been aligned (using Linux tools)

os.system("samtools stats aln.sorted.sam")

