#!/bin/bash

#check if necessary arguments are introduced
if [ $# -ne 2 ]; then
 echo "The program needs to arguments, genome file and reads file (in this order)";
 exit;
fi
GENOME=$1
FILE=$2

# Make temp directory for files produced
if [ ! -d assign4_temp ]; then
  mkdir assign4_temp;
fi

FIRSTCHAR=$(cat $FILE | head -c 1);
# If first character is > run fasta script
if [ $FIRSTCHAR = ">" ]; then
    printf "\nInput sequence identified as fasta type\n"
    printf "\nReading input sequence, counting number of 3-mers ('mer-counts.txt')...\n";
## 1. Counts the abundance of each possible 3-mers (histogram)
    grep -v "^>" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 3 > trimers.1.txt;
    grep -v "^>" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 1 | tail -n +2 | tr -d '\n' | fold -w 3 > trimers.2.txt;
    grep -v "^>" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 1 | tail -n +3 | tr -d '\n' | fold -w 3 > trimers.3.txt;
    
## 2. Splits the input in FASTA/FASTQ into files of only 1 sequence each
printf "\nCreating seperate fasta files...\n";
# Split to seperate fasta files
#head -n 20000 $FILE | awk '/>/{x="assign4_temp/fasta_seq"++i".fasta";}{print > x;}';
cat  $FILE | awk '/>/{x="assign4_temp/fasta_seq"++i".fasta";}{print > x;}';

##3. Hard-trims (from the right) all sequences from all files 20nt
# loop to make all nucleotides uppercase
for file in assign4_temp/*.fasta; do
    head=$(cat $file | head -n 1); # Leave head as is
    rest=$(cat $file | tail -n +2 $file | tr '[:lower:]' '[:upper:]'); # The rest transform to uppercase
    rest_no_whitespace="$(echo -e "${rest}" | tr -d '[:space:]')";
    rest_cut=${rest_no_whitespace::-20}; # Trim last 20 bases
    logwrite="$head\n$rest_cut";
    echo -e $logwrite > $file;
done
fi

# If first character is @ run fastq script
if [ $FIRSTCHAR = "@" ]; then
    printf "\nInput sequence identified as fastq type\n"
    printf "\nReading input sequence, counting number of 3-mers ('mer-counts.txt')...\n";
## 1. Counts the abundance of each possible 3-mers (histogram)
    grep "[ACTGN]" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 3 > trimers.1.txt;
    grep "[ACTGN]" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 1 | tail -n +2 | tr -d '\n' | fold -w 3 > trimers.2.txt;
    grep "[ACTGN]" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 1 | tail -n +3 | tr -d '\n' | fold -w 3 > trimers.3.txt;


## 2. Splits the input in FASTA/FASTQ into files of only 1 sequence each
    printf "\nCreating seperate fastq files...\n";
# Split to seperate fasta files
    split -l 4 -d --additional-suffix=.fastq $2 assign4_temp/fastq_seq;

##3. Hard-trims (from the right) all sequences from all files 20nt
# loop to make all nucleotides uppercase
    for file in assign4_temp/*.fastq; do
        head=$(cat $file | sed -n 1p); # Leave head as is
        seq=$(cat $file | sed -n 2p $file | tr '[:lower:]' '[:upper:]'); # The rest transform to uppercase
        qual=$(cat $file | sed -n 4p $file); # The rest transform to uppercase
        seq_cut=$(echo $seq | head -c -20); # Trim last 20 bases
        qual_cut=$(echo $qual | head -c -20); # Trim last 20 bases
        echo $head > $file ;
        echo $seq_cut >> $file ;
        echo '+' >> $file ;
        echo $qual_cut >> $file ;
    done
fi



# Read stats file, only keep the nucleotides (some files had 3-mers like AMW or WGA for some reason) common for fasta and fastq
cat trimers.1.txt trimers.2.txt trimers.3.txt | sort | uniq -c > stats.txt;
    rm trimers.1.txt;
    rm trimers.2.txt;
    rm trimers.3.txt;
    cat stats.txt| sort | awk -F ' ' 'length($2)== 3{print$1,$2}' > stats_filt.txt;
    rm stats.txt;

echo > mer-counts.txt;
cat stats_filt.txt | while read line 
do
  line=$(echo $line | tr '[:lower:]' '[:upper:]' )
  IFS=$' '; # space is set as delimiter
  read -ra ADDR <<< "$line"; # str is read into an array as tokens separated by IFS
  mer="${ADDR[1]}"; # get mer
  valid='[ACTG][ACTG][ACTG]'
  if [[ $mer =~ $valid ]]; then # If mer only contains nucleotide bases, append to stats_filtered
    echo $line >> mer-counts.txt;
  fi
done

rm stats_filt.txt;


# Keep only sam files that produced alignment
printf "\nTransforming SAM files...\n";
for file in assign4_temp/*.sam; do
  heads="$(cat $file | grep "^@")"; # read file
  line="$(cat $file | grep -v "^@")";
  IFS=$'\t'; # tab is set as delimiter
  read -ra ADDR <<< "$line"; # str is read into an array as tokens separated by IFS
  chr="${ADDR[2]}"; # get chromosome id
  start="${ADDR[3]}"; # get start position in sequence
  rm $file; # remove old sam file
  # if chromosome id and start position are different than '*' keep the sam file
  if [ "$chr" != "*" ]; then
    if [ "$start" != "*" ]; then
      echo -e "$heads\n$line" > $file;
    fi
  fi
done

##4. Using the genome mapping tool BWA and the reference genome of the Scaromice Cerevisiae (any strain will do)
## aligns each of the files producing its corresponding SAM file.
printf "\nIndexing the reference sequence and creating alignment files...\n";
# index genome
bwa index -p GENIDX -a bwtsw $GENOME 2> /dev/null; 
# for each fastq create sam allignment file
for file in assign4_temp/*.fast*; do bwa mem GENIDX $file > $file.sam 2> /dev/null; rm $file ;done 


#5. Merge all SAM files ignoring headers (using Linux tools)
#6. Sorts the SAM file by chromosome and position
# Sort and merge sam files
printf "\nSorting and merging SAM files...\n"
for file in assign4_temp/*.sam; do
  samtools sort -O SAM -o $file $file;
done;
samtools merge -O sam -f merge.sam assign4_temp/*.sam > merge.sam;
for file in assign4_temp/*.sam; do rm $file; done; # remove sam files from assign4_temp dir
cat merge.sam | grep -v "^@" > sorted-merged.sam; # remove header lines


printf "\nRemoving files produced by the program...\n\n";
# rm temp directory
if [ ! "$(ls -A "assign4_temp")" ];
  then
    rmdir assign4_temp;
fi
rm GENIDX.*;

printf "\nDone\n";

#7. Compute how many reads have been aligned (using Linux tools)
# Count number of mapped reads
printf "\nNumber of reads mapped:\n";
samtools view merge.sam | wc -l;
rm merge.sam;

printf "\nThe sorted merged SAM file can be found in the direcory under the name:\nsorted-merged.sam:\n";



