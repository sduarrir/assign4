import sys
##FUNCTIONS

# Function to open input file and check if correct format
def open_input(filename): #we provide the file
    try:
        inputfile = open(filename, "rt")
    except FileNotFoundError:
        print(filename, "Not found. check for any typos!")
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

#trims and also counts trinucleotides
def proc_seq(seq, kdict, rtrim):
    for i in range(len(seq) - 2):
        kmer = sequence[i:i + 3]
        if kmer not in kdict:
            kdict[kmer] = 1
        else:
            kdict[kmer] += 1
    return seq[:-rtrim], dict

#check if arguments were checked
if len(sys.argv) != 3:
    print("The program need the following arguments: file where geneome is stored and file where the the reads are stored\nFormat must be fasta or fastq\n")
    exit(1)

#check if fasta or fastaq
genomefile, _ = open_input(sis.argv[1])
genomefile.close()
readsfile, format = open_input(sis.argv[1])
statfile= open(sys.argv[1] + ".stats","wt")
kdict = {}
cut = 20
i = 0

#ouput files should be on the temp folder!!!
while True:
    tag = readsfile.readline()
    if not tag: break
    #read seq
    outf = open('read' + i + '.' + format, 'wt')
    seq = readsfile.readline.strip()
    seq, kdict = proc_seq(seq, kdict, cut)
    outf.write("%s%s\n" % (tag, seq))
    if format == 'fastq':
        plus = readsfile.readline()
        qual = readsfile.readline.strip()
        qual = qual[:-cut]
        outf.write("%s%s\n" % (plus, qual))
    outf.close()


input_file.close()

# Open stats file
output_file= open(sys.argv[1] + ".stats","wt")
output_file.write(str(occ))
output_file.close()
