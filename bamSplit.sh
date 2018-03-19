#!/bin/bash
# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 2.0
# 2015-07-13
# 2015-09-02
#
# bamSplit.sh

########################################################
#                                                      #
#                                                      #
#   This script split BAM files by chromosome          #
#                                                      #
#                                                      #
########################################################

menu="
bamSplit.sh version 2.0
--------------------------------
This script split BAM files by chromosome.

Usage   : ./bamSplit.sh inFile.bam outDirectoryName
Example : ./bamSplit.sh IP_exp1.bam experiment1

This script will generate an output similar to : experiment1/IP_exp1.REF_*.bam
Where \"*\" represents each chromosome.
"

# Display color in shell
ERROR="`tput setaf 1`[ERROR]`tput sgr0`"

# Display help
if [[ $1 == -h ]] || [[ $1 == --help ]]
then
    echo "$menu"
    exit
fi

# Validation of user's input
if ! [ $# -eq 2 ]
then
    echo ""
    echo " $ERROR : Incorrect amount of arguments : requires 2, $# provided."
    echo "$menu"
    exit
fi

# Validation of input files (must exists)
if [ ! -f $1 -o -h $1 ]
then
    echo " $ERROR : File $1 not found !"
    exit
fi

# Assign input parameters
bamFile=`readlink -f $1`
outDir=$2

# Validation of input extension
if [[ ! ${bamFile##*.} == bam ]]
then
    echo " $ERROR : File extension not accepted, BAM files only !"
    exit
fi

# Validation of output directory (must not exists)
if [ -d $outDir ]
then
    echo " $ERROR : Directory $outDir exists !"
    read -p "Do you want to replace existing directory ? [yes|no] : " choice
    
    until [[ $choice == yes ]] || [[ $choice == no ]]
    do
        echo " $ERROR : Please answer yes or no !"
        read -p "Do you want to replace existing directory ? [yes|no] : " choice
    done
    
    case $choice in

        yes)
            echo "Directory $outDir will be overwritten !"
            rm -r $outDir
            ;;
            
        no)
            echo "Please, change your output directory name."
            echo "Exiting program..."
            exit
            ;;
        *)
            echo "`tput setaf 1`[ERROR101]`tput sgr0` : An internal error occured, \
            please report it to marianne.sabourin-felix.1@ulaval.ca"
            exit
            ;;
    esac
fi

StartTime=$(date +%s)

# Create new directory
mkdir $outDir

# Split BAM files
bamtools split -in $bamFile -reference

# Move split files to output directory
path="${bamFile%%.*}"
mv $path.REF_*.bam $outDir

# Replace .REF by -REF
for file in $outDir/*
do
    mv $file ${file//\.REF/-REF}
done


EndTime=$(date +%s)
ElapsedTime=$(($EndTime - $StartTime))

echo "Splitting done for file `basename $bamFile` in $(($ElapsedTime / 3600 )) \
hours $((($ElapsedTime % 3600) / 60)) minutes $(($ElapsedTime % 60)) seconds !"

