#!/bin/bash
##ASSIGNMENT 4

#check if necessary arguments are introduced
if [ $# -ne 2 ]; then
 echo "The program needs to arguments, genome file and reads file (in this order)";
 exit;
fi


## 1. Counts the abundance of each possible 3-mers (histogram)
# Make temp directory for files produced
if [ ! -d assign4_temp ]; then
  mkdir assign4_temp;
fi

#save arguments
FILE=$2; #FILE WITH THE READS IS THE SECOND ARGUMENT
GENOME=$1
# Count 3-mers in file
grep -v "^>" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 3 > trimers.1.txt; #orf1
grep -v "^>" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 1 | tail -n +2 | tr -d '\n' | fold -w 3 > trimers.2.txt; #orf2
grep -v "^>" $FILE | tr '[:lower:]' '[:upper:]' | fold -w 1 | tail -n +3 | tr -d '\n' | fold -w 3 > trimers.3.txt; #orf3
cat trimers.1.txt trimers.2.txt trimers.3.txt | sort | uniq -c > stats.txt; #sorted trimers
rm trimers.1.txt;
rm trimers.2.txt;
rm trimers.3.txt;

## remove those that are not trimers len < 3 (max number of lines of stats should be 64)

cat stats.txt | awk -F ' ' 'length($2)== 3{print$1,$2}' > stats_filt.txt # now need to save it in a file!!

## 2. Splits the input in FASTA/FASTQ into files of only 1 sequence each
printf "\nCreating seperate fasta files...\n";
# Split to seperate fasta files 
csplit sequence.fasta '/>/' {*} --prefix='./assign4_temp/seq' --suffix-format='%d.fasta' --elide-empty-files -s #silent
#this could lead into problems if in the header there is a > symbol



##3. Hard-trims (from the right) all sequences from all files 20nt
# loop to make all Nucleotides uppercase #changed the way of cutting (we avoid trimmind after id newline)
for file in assign4_temp/*; do
    head=$(cat $file | head -n 1); # Leave head as is
    rest=$(cat $file | tail -n +2 $file | tr '[:lower:]' '[:upper:]'); # The rest transform to uppercase
    rest_no_whitespace="$(echo -e "${rest}" | tr -d '[:space:]')"
    rest_cutted=${rest_no_whitespace::-20} #TRIM LAST 20!
    echo $rest_cutted >> $file; # Append to file
done

##4. Using the genome mapping tool BWA and the reference genome of the Scaromice Cerevisiae (any strain will do), aligns each of the files producing its corresponding SAM file.
printf "\nIndexing Saccharomyces (S288C) genome and creating alignment files...\n";
 # index sacch genome sacch as db index(p) indicate we will be woroking with a whole genome(-a bwtsw)
bwa index -p sacch -a bwtsw $GENOME 2> /dev/null;

# for each fasta create sam allignment file
#for file in assign4_temp/*.fasta; do bwa bwasw sacch $file > $file.sam 2> /dev/null ;rm $file;done
for file in assign4_temp/*.fasta; do bwa bwasw sacch $file > $file.sam 2> /dev/null ;done

#5. Merge all SAM files ignoring headers (using Linux tools)
printf "\nRemoving SAM headers and creating a merged SAM file (merged-file.sam)...\n";
# creates merged sam file from all the separate sam files in assign4_temp

echo -n > merged-file.sam; # open new file "merged-file.sam" !!!! THIS WILL DELETE PREVIOUS FILES!!!!
for file in assign4_temp/*.sam; do
  str="$(cat $file)"; # read file
  str="$(grep -v '^@' <<<$str)"; # remove header lines
  IFS=$'\t'; # space is set as delimiter
  read -ra ADDR <<< "$str"; # str is read into an array as tokens separated by IFS
  cigar="${ADDR[5]}"; # Look for the cigar str in SAM
  
  if [ "$cigar" != "*" ]; then # If cigar is not empty
    echo "$str" >> merged-file.sam; # append the sam file to merged-sam file
  fi
  #rm $file; # remove the file
done;
#SAM HAS A TOOL THAT FILERS NON ALIGNED SEQUENCES (so we could do this after the merge!##


#6. Sorts the SAM file by chromosome and position
samtools sort merged-file.sam > sorted-merged-file.sam #still NOT WORKS :((8
#3rd colum indicates th cromosome (i think they're indexed) 4th starting position

#7. Compute how many reads have been aligned (using Linux tools)
wc -l sorted-merged-file.sam

##THIS IS MISSING


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


