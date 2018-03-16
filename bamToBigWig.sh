#!/bin/bash

# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 1.0
# 2018-03-15
#
# bamToBigWig.sh

# Custom display in shell
red=`tput setaf 1`
bold=`tput bold`
normal=`tput sgr0`
ERROR="$red[ERROR]$normal"

# Help
menu="
$bold--------------------------------------------------------
 bamToBigWig.sh version 1.0
--------------------------------------------------------$normal
 This script converts bam files into bigWig

   Usage   : ./bamToBigWig.sh file.bam /path/to/bedGraphToBigWig

$bold [Dependecies] $normal: - Bedtools
                 - bedGraphToBigWig (UCSC)
$bold--------------------------------------------------------$normal
"

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
if ! [ $# -eq 2 ]
then
    echo ""
    echo "  $ERROR Incorrect amount of arguments : requires 2, $# provided."
    echo "$menu"
    exit
fi

# Assignation of input parameter
bamFile=$1
pathToBedGraphToBigWig=$2

# Validation of input file
if ! [ -f $bamFile -o -h $bamFile ]
then
    echo " $ERROR : File" `basename $bamFile` "not found !"
    exit
fi

# Validation of input extension
if [[ $bamFile != *.bam ]]
then
	echo " $ERROR : File" `basename $bamFile` "must be in BAM format !"
    exit
fi

# Validation if path to folder exists
if ! [ -d $pathToBedGraphToBigWig -o -L $pathToBedGraphToBigWig ]
then
    echo " $ERROR : Folder $pathToBedGraphToBigWig not found !"
    exit
fi

# Validation of bedGraphToBigWig path
if ! [[ `find $pathToBedGraphToBigWig/. -maxdepth 1 -type f -name "bedGraphToBigWig" 2>/dev/null` ]]
then
	echo " $ERROR : The folder $pathToBedGraphToBigWig must contain the bedGraphToBigWig program !"
    exit
fi


##################################
#                                #
#########     M A I N    #########
#                                #
##################################

# Creation of the chromInfo file
samtools idxstats $bamFile | cut -f 1,2 > chromInfo.txt

# Creation of output filenames
outputBdg=`basename $bamFile .bam`.bedGraph
outputBw=`basename $bamFile .bam`.bw

# Conversion BAM -> BedGraph
bedtools genomecov -bg -ibam $bamFile > $outputBdg

# Conversion BedGraph -> BigWig
$pathToBedGraphToBigWig/bedGraphToBigWig $outputBdg chromInfo.txt $outputBw

# Remove intermediate file
rm $outputBdg
rm chromInfo.txt

