#!/bin/bash
##ASSIGNMENT 4


3. Hard-trims (from the right) all sequences from all files 20nt
4. Using the genome mapping tool BWA and the reference
genome of the Scaromice Cerevisiae (any strain will do), aligns
each of the files producing its corresponding SAM file.
5. Merge all SAM files ignoring headers (using Linux tools)
6. Sorts the SAM file by chromosome and position
7. Compute how many reads have been aligned (using Linux
tools)
21


## 1. Counts the abundance of each possible 3-mers (histogram)
# Make temp directory for files produced
if [ ! -d assign4_temp ]; then
  mkdir assign4_temp;
fi

# Count 3-mers in file
FILE=sequence.fasta;
grep -v "^>" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 3 > trimers.1.txt; #orf1
grep -v "^>" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 1 | tail -n +2 | tr -d '\n' | fold -w 3 > trimers.2.txt; #orf2
grep -v "^>" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 1 | tail -n +3 | tr -d '\n' | fold -w 3 > trimers.3.txt; #orf3
cat trimers.1.txt trimers.2.txt trimers.3.txt | sort | uniq -c > stats.txt; #sorted trimers
rm trimers.1.txt;
rm trimers.2.txt;
rm trimers.3.txt;

## remove those that are not trimers len < 3 (max number of lines of stats should be 64)



## 2. Splits the input in FASTA/FASTQ into files of only 1 sequence each
printf "\nCreating seperate fasta files...\n";
# Split to seperate fasta files 
csplit ../sequence.fasta '/>/' {*} --prefix='/assign4_temp/seq' --suffix-format='%d.fasta' --elide-empty-files 
#this could lead into problems if in the header theres a > symbol

#head -n 1600 sequence.fasta | awk '/>/{x="assign4_temp/fasta_seq"++i".fasta";}{print > x;}';
#that just created a file dit not split as it suposed to

# loop to make all Nucleotides uppercase #changed the way of cutting (we avoid trimmind after id newline)
for file in assign4_temp/*; do
    head=$(cat $file | head -n 1); # Leave head as is
    rest=$(cat $file | tail -n +2 $file | tr '[:lower:]' '[:upper:]'); # The rest transform to uppercase
    rest_no_whitespace="$(echo -e "${rest}" | tr -d '[:space:]')"
    rest_cutted=${rest_no_whitespace::-20} #TRIM LAST 20!
    echo $rest_cutted >> $file; # Append to file
done


printf "\nIndexing Saccharomyces (S288C) genome and creating alignment files...\n";
 # index sacch genome
bwa index -p sacch -a bwtsw sacch.fasta 2> /dev/null;

# for each fasta create sam allignment file
for file in assign4_temp/*.fasta; do bwa bwasw sacch $file > $file.sam 2> /dev/null ;rm $file;done 

printf "\nRemoving SAM headers and creating a merged SAM file (merged-file.sam)...\n";
# creates merged sam file from all the separate sam files in assign4_temp
for file in assign4_temp/*.sam; do
  str="$(cat $file)"; # read file
  str="$(grep -v '^@' <<<$str)"; # remove header lines
  IFS=$'\t'; # space is set as delimiter
  read -ra ADDR <<< "$str"; # str is read into an array as tokens separated by IFS
  cigar="${ADDR[5]}"; # Look for the cigar str in SAM
  echo -n > merged-file.sam; # open new file "merged-file.sam"
  if [ "$cigar" != "*" ]; then # If cigar is not empty
    echo "$str" >> merged-file.sam; # append the sam file to merged-sam file
  fi
  rm $file; # remove the file
done;

##################### STILL NEED TO SORT THE MERGED SAM FILE #############################

printf "\nRemoving files from direcory produced by the program...\n\n";
# rm temp directory
if [ "$(ls -A "assign4_temp")" ]
  then
     pass;
  else
    rmdir assign4_temp;
fi

# remove indexed sacch files produced by bwa
for file in *sacch*; do
  if [ "$file" != "sacch.fasta" ]
  then
    rm $file;
  fi
done;
printf "\nDone\n";
printf "\nFirst 10 lines from the sorted-merged-file.sam:\n\n";
head sorted-merged-file.sam;




