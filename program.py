"""
            ---     TASK    ---
Produce two scripts (Bash/Python) that, given an input FASTA/FASTQ file, implement the following pipeline:
"""
import sys
import os
#from Bio import SeqIO #would allow us to process sequences
##FUNCTIONS


# Function to open input file and check if correct format
def open_input(filename): #we provide the file
    try:
        inputfile = open(filename, "rt")
    except FileNotFoundError:
        print(filename, "not found. check for any typos!")
        exit(1)
    line = inputfile.readline()
    if line[0] == "@":
        fileformat = "fastq"
        print("%s format is fastq" % filename)
    elif line[0] == ">":
        fileformat = "fasta"
        print("%s format is fasta" % filename)
    else:
        print(filename, "format is not fastaq or fasta.") #if format is not fasta or fastaq it will print an error and exit
        exit(1)
    inputfile.seek(0)
    return (inputfile, fileformat) #returns the file and its format (fasta o fastaq)

#trims and also counts trinucleotides
def proc_seq(seq, kdict, rtrim):
    for i in range(len(seq) - 2):
        kmer = seq[i:i + 3]
        if kmer not in kdict:
            kdict[kmer] = 1
        else:
            kdict[kmer] += 1
    return seq[:-rtrim], kdict

def rm_files(files): #given a file / path removes what is indicated
    clean = 'rm ' + files
    os.system(clean)



#open the files (and check if format is correct
if len(sys.argv) !=3:
    print("The program need the following arguments:\n\t--input [INPUTFILE]\n\t--ouput [FILE]  If file exists, will be overwrited!!!\n\t--operation [OPERATION] (what do you want the program to print\nSome operations will need conditions.\n\trc (reverse complement) does not need conditions\n\ttrim needs: --trim-right [POSITION] and --trim-left [POSITION]\n\tadaptor-removal needs: --adaptor [SEQUENCE]")
    exit(1)
genome = sys.argv[1]
genomefile, _ = open_input(genome)
genomefile.close()
reads = sys.argv[2]
readsfile, format = open_input(reads)

kdict = {}
cut = 20
i = 0 #count for the files

##1.Counts the abundance of each possible 3-mers (histogram)
##2.Splits the input in FASTA/FASTQ into files of only 1 sequence each
##3.Hard-trims (from the right) all sequences from all files 20nt
path = "assignment4_temp/"
clean = path + '*' #to clean path
try:
    os.mkdir(path)
except:
    print ('Directory already exists. Removing files')
    rm_files(clean)

print("\nReading input sequence, counting number of 3-mers ('mer-counts.txt')...\nCreating seperate fastx files...\n")
while True:
    tag = readsfile.readline()
    if not tag: break
    #read seq
    filename = path + 'read' + str(i) + '.' + format
    outf = open(filename, 'wt')
    seq = readsfile.readline().strip()
    seq, kdict = proc_seq(seq.upper(), kdict, cut) #change  sequence top uppercase
    outf.write("%s%s\n" % (tag, seq))
    if format == 'fastq':
        plus = readsfile.readline()
        qual = readsfile.readline().strip()
        qual = qual[:-cut]
        outf.write("%s%s\n" % (plus, qual))
    outf.close()
    i += 1
readsfile.close()

statfile = open(sys.argv[1] + ".stats","wt")

for key in sorted(kdict):
    statfile.write("%s : %d" % (key, kdict[key]))
#statfile.write(kdict)


##4.Using the genome mapping tool BWA and the reference genome of the ScaromiceCerevisiae(any strain will do), aligns each of the files producing its corresponding SAM file.
gencode = 'SACCHGEN'
command = 'bwa index -p ' + gencode + ' -a bwtsw ' + genome + ' 2> /dev/null'
os.system(command)
print("\nIndexing the reference sequence and creating alignment files...\n")
for file in os.listdir(path):
    file = path + file
    file2 = file + '.sam'
    f=open(file2, 'wt')
    f.write('@HD\tVN:1.0\tSO:unsorted\n')
    f.close()
    command = 'bwa mem -a ' + gencode + ' ' + file + ' 1>> ' + file2 + ' 2> /dev/null'
    os.system(command)
    os.remove(file)

##5.Merge all SAM files ignoring headers (using Linux tools)
merged = 'merged.sam'
command = 'samtools merge -O sam -f ' + merged + ' ' + path + '*.sam'
os.system(command)
#rm_files(clean)
#os.system(command)
#TILL HERE IT WORKS (MERGED SAM FILE)
#rm_files(file)

##6.Sorts the SAM file by chromosome and position
print("\nTransforming SAM file...\n")
fmerged = open(merged, 'rt')
filt = 'filt_merge.sam'
ffilt = open(filt, 'wt')
count = 0
while True: #keep the c
    read = fmerged.readline()
    if not read : break
    if read[0] == '@':
        #if read[:3] != '@HD':
            #continue
        ffilt.write(read)
    else:
        readlist = read.split()
        if readlist[1] != '*' and readlist[2] != '*':
            count += 1
            read = '\t'.join(readlist) + '\n'
            ffilt.write(read)

fmerged.close()
rm_files(merged)
ffilt.close()

command = 'samtools sort -O SAM -o '+ filt+ ' ' + filt
os.system(command)
f = open('final.sam', 'wt')
ffilt=open(filt, 'rt')
#clean output file
while True:
    read = ffilt.readline()
    if not read : break
    if read[0] != '@':
        f.write(read)
ffilt.close()
rm_files(filt)
f.close()

##7.Compute how many reads have been aligned (using Linux tools)

print('%s reads aligned' % count)

