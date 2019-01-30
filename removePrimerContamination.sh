#!/bin/bash

# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 1.0
# 2016-02-10
#
# removePrimerContamination.sh

# Help
menu="
removePrimerContamination.sh version 1.0
-----------------------------------
This script removes primer contamination from BAM files.

Usage   : ./removePrimerContamination.sh processedFile1.bam processedFile2.bam processedFile3.bam processedFile4.bam outputFolder

 processedFile*.bam -> Sorted alignment file (Input)

This script will generate output similar to : outputFolder/processedFile*_corrected.bam
"

# Display color in shell
ERROR="`tput setaf 1`[ERROR]`tput sgr0`"

# Display help
if [[ $1 == -h ]] || [[ $1 == --help ]]
then
    echo "$menu"
    exit
fi


##################################
#                                #
# Validation of user's input     #
#                                #
##################################

# Validation of input parameters
if ! [ $# -eq 5 ]
then
    echo ""
    echo "  $ERROR Incorrect amount of arguments : requires 5, $# provided."
    echo "$menu"
    exit
fi

# Assign input parameters
file1=$1
file2=$2
file3=$3
file4=$4
outputFolder=$5

# Validation of input file
for file in $file1 $file2 $file3 $file4
do
    if ! [ -f $file -o -h $file ]
    then
        echo " $ERROR : File" `basename $file` "not found !"
        exit
    fi
done

# Validation of input extension
for file in $file1 $file2 $file3 $file4
do
    if [[ ! ${file##*.} == bam ]]
    then
        echo " $ERROR : File extension not accepted, BAM files only !"
        exit
    fi
done

# Validation of output directory (must not exists)
if [ -d $outputFolder ]
then
    echo " $ERROR : Directory $outputFolder exists !"
    read -p "Do you want to replace existing directory ? [yes|no] : " choice
    
    until [[ $choice == yes ]] || [[ $choice == no ]]
    do
        echo " $ERROR : Please answer yes or no !"
        read -p "Do you want to replace existing directory ? [yes|no] : " choice
    done
    
    case $choice in

        yes)
            echo "Directory $outputFolder will be overwritten !"
            rm -r $outputFolder
            ;;
            
        no)
            echo "Please, change your output name."
            echo "Exiting program..."
            exit
            ;;
        *)
            echo "`tput setaf 1`[ERROR101]`tput sgr0` : An internal error occured, please report it to marianne.sabourin-felix.1@ulaval.ca"
            exit
            ;;
    esac
fi


##################################
#                                #
# Remove primers contamination   #
#                                #
##################################

StartTime=$(date +%s)

# Creation of output directory
mkdir $outputFolder

primerForward=GTGAGCGATCTGTATTGGTTTGTATGGTTGATCGAGACCATTGTCG
primerReverse=GGGCAAGCAGGCGGGAGCGTCTCGGAGATGGTGTCGTGTTTAAGGA

for file in $file1 $file2 $file3 $file4
do
    echo "Processing file" `basename $file`"..."
    
    # Extract filename
    prefix=`basename $file .bam`
    
    # Conversion bam -> sam
    echo "Conversion bam -> sam..."
    samtools view -h $file > $outputFolder/$prefix.sam
    
    # Detection of forward primer
    nbForward=`grep "$primerForward" $outputFolder/$prefix.sam | wc -l`
    echo "The forward primer was found $nbForward times in file $file"
    
    # Detection of reverse primer
    nbReverse=`grep "$primerReverse" $outputFolder/$prefix.sam | wc -l`
    echo "The reverse primer was found $nbReverse times in file $file"
    
    # Remove primer sequences
    echo "Removing primer sequences..."
    grep -vE "$primerForward|$primerReverse" $outputFolder/$prefix.sam > $outputFolder/$prefix"_corrected.sam"
    
    # Remove not corrected sam file
    rm $outputFolder/$prefix.sam
    
    # Conversion sam -> bam
    echo "Conversion sam -> bam..."
    samtools view -bS $outputFolder/$prefix"_corrected.sam" > $outputFolder/$prefix"_corrected.bam"
    
    # Remove corrected sam file
    rm $outputFolder/$prefix"_corrected.sam"
    
    # Index of final bam file
    echo "Indexing final bam file..."
    samtools index $outputFolder/$prefix"_corrected.bam"
    
    echo "Final file" $prefix"_corrected.bam created !"

done

EndTime=$(date +%s)
ElapsedTime=$(($EndTime - $StartTime))

echo "Primer removed for experiment $outputFolder in $(($ElapsedTime / 3600 )) \
hours $((($ElapsedTime % 3600) / 60)) minutes $(($ElapsedTime % 60)) seconds !"




