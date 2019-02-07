#!/bin/bash

# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 1.0
# 2017-03-23
#
# GCcontentOfPeaksCalledByMACS2.sh

###############################################
#                                             #
#  This script allows to extract the GC% of   #
#  peaks called by MACS2.                     #
#                                             #
###############################################


# Help
menu="
GCcontentOfPeaksCalledByMACS2.sh version 1.0
------------------------------------------------------
This script allows to extract the GC% of peaks called by MACS2.

Usage : ./GCcontentOfPeaksCalledByMACS2.sh narrowPeak chrFolder outputName

narrowPeak -> Peak files called with MACS2 (narrowPeak format)
chrFolder  -> Folder that contain the chromosomes files in fasta with the sequence on one line without header *
outputName -> Output name of the GC% file
 
* To put the sequence on one line do : sed '1d' file.fa | tr -d '\n' output.fa

This script will generate output similar to : outputName_GC.tsv
that contain : chr, start, end, peakName, sequence, GC%
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
#   Validation of user's input   #
#                                #
##################################

# Validation of input parameters
if ! [ $# -eq 3 ]
then
    echo ""
    echo "  $ERROR Incorrect amount of arguments : requires 3, $# provided."
    echo "$menu"
    exit
fi

# Assignation of input parameter
narrowPeak=$1
chrFolder=$2
outputName=$3

# Validation of input file
if ! [ -f $narrowPeak -o -h $narrowPeak ]
then
    echo " $ERROR : File" `basename $narrowPeak` "not found !"
    exit
fi

# Validation of input extension
if [[ $narrowPeak != *.narrowPeak ]]
then
	echo " $ERROR : File" `basename $narrowPeak` "must be in narrowPeak format !"
    exit
fi

# Validation of number of columns (must be 10)
nbOfCol=`awk '{print NF; exit}' $narrowPeak`

if [[ $nbOfCol != 10 ]]
then
	echo " $ERROR : File" `basename $narrowPeak` "must contain 10 columns !"
	exit
fi

# Validation of chromosome folder
if ! [[ `find $chrFolder/. -maxdepth 1 -type f -name *.fa 2>/dev/null` ]]
then
	echo " $ERROR : The folder $chrFolder must contain *.fa files !"
    exit
fi

# Creation of output filename
outputFile=$outputName".tsv"

# Validation of output file (must not exists)
if [ -f $outputFile -o -h $outputFile ]
then
    echo " $ERROR : File" `basename $outputFile` "exists !"
    read -p "Do you want to replace existing file ? [yes|no] : " choice
    
    until [[ $choice == yes ]] || [[ $choice == no ]]
    do
        echo " $ERROR : Please answer yes or no !"
        read -p "Do you want to replace existing file ? [yes|no] : " choice
    done
    
    case $choice in

        yes)
            echo "File $outputFile will be overwritten !"
            rm -r $outputFile
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
#########     M A I N    #########
#                                #
##################################

StartTime=$(date +%s)

# Read the narrowPeak file line by line and extract the columns into parameters
while read -r chr start end peakName score strand signalValue pValue qValue pointSource
do
	# Extract the sequence of each peak
	newStart=$((start + 1))
	len=$((end - start))
	sequence=`awk -v a=$newStart -v b=$len '{print substr($0, a, b)}' $chrFolder/$chr.fa`

	# Compute the GC% for each peak from its sequence
    sumCount=0
	for base in `echo $sequence | fold -1`
	do

		case $base in
    
        	C|c|G|g|S|s)
            	# S -> C or G
            	count=100
            	;;
        
        	B|b|V|v)
            	# B -> C or G or T
            	# V -> A or C or G
            	count=66.666
            	;;
        
        	K|k|M|m|R|r|Y|y)
            	# K -> G or T
            	# M -> A or C
            	# R -> A or G
            	# Y -> C or T
            	count=50
            	;;
        
        	D|d|H|h)
            	# D -> A or G or T
            	# H -> A or C or T
            	count=33.333
            	;;
        
        	A|a|T|t|W|w|N|n)
            	# W -> A or T
            	# N -> A or C or G or T
            	count=0
            	;; 
        
        	*)
            	echo "`tput setaf 1`[ERROR102]`tput sgr0` : An internal error occured, please report it to marianne.sabourin-felix.1@ulaval.ca"
            	exit
            	;;
        esac

        sumCount=$((sumCount+count))

	done

	gc=$((sumCount/len))

	echo -e "$chr\t$start\t$end\t$peakName\t$sequence\t$gc" >> $outputFile

done < $narrowPeak

mean=`awk '{ total += $6 } END { print total/NR }' $outputFile`

echo "The mean of all the peaks is : $mean %"

EndTime=$(date +%s)
ElapsedTime=$(($EndTime - $StartTime))

echo "--------------------------------------------------------------------------------"

echo "Sequence extraction and GC% calculation done in $(($ElapsedTime / 3600 )) \
hours $((($ElapsedTime % 3600) / 60)) minutes $(($ElapsedTime % 60)) seconds !"


