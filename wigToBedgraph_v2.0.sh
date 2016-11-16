#!/bin/bash
# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 2.0
# 2015-07-08
# 2015-09-02
#
# wigToBedgraph.sh

########################################################
#                                                      #
#                                                      #
#   This script convert wig files to bedgraph files    #
#                                                      #
#                                                      #
########################################################

menu="
wigToBedgraph.sh version 2.0
----------------------------------
This script convert wig files to bedgraph files.

Usage   : ./wigToBedgraph.sh inFile.wig outFileName
Example : ./wigToBedgraph.sh MmrDNA35.wig experiment1

This script will generate output similar to : experiment1.bedgraph
  --> If the output file exists, it will be overwritten !
"

# Display color in shell
ERROR="`tput setaf 1`[ERROR]`tput sgr0`"

# Display help
if [[ $1 == -h ]] || [[ $1 == --help ]]
then
    echo "$menu"
    exit
fi

# Validation of input user
if ! [ $# -eq 2 ]
then
    echo ""
    echo " $ERROR : Incorrect amount of arguments : requires 2, $# provided."
    echo "$menu"
    exit
fi

# Validation of input files
if [ ! -f $1 -o -h $1 ]
then
    echo " $ERROR : File $1 not found !"
    exit
fi

# Assign input parameters
wigFile=$1
outName=$2
outputFile=$outName.bedgraph

# Validation of input extension
if [[ ! ${wigFile##*.} == wig ]]
then
    echo " $ERROR : File extension not accepted, wig files only !"
    exit
fi

# fixedStep not supported yet...
if [ "`grep 'fixedStep' $wigFile`" ]
then
    echo " $ERROR : fixedStep wig format not supported yet !"
    echo "Please report this error to marianne.sabourin-felix.1@ulaval.ca"
    exit
fi

StartTime=$(date +%s)

# Print header into file (this step will overwrite an existing file with the same name)
echo -e "track type=bedgraph name=\"$outName mappability\" description=\"User Supplied Track\" visibility=full color=0,0,255 altColor=255,0,0 priority=20 autoScale=on alwaysZero=off gridDefault=off maxHeightPixels=128:128:11 graphType=bar viewLimits=0.0:25.0 yLineMark=0.0 yLineOnOff=off windowingFunction=maximum smoothingWindow=off" > $outputFile

step=1
# sed '1d' $wigFile | while read -r line
while read -r line
do
    # Skip trackline
    [ "`echo $line | grep 'track'`" ] && continue
    
    ### fixedStep - Not supported yet... ###
    # format : fixedStep chrom=chr3 start=400601 step=100 span=5
    # format : dataValue
    
    ### VariableStep ###
    # Get chromosome name
    if [ "`echo $line | grep 'chrom'`" ]
    then
        chromName=`echo $line | grep -o 'chrom.*' | cut -f1 -d' ' | cut -f2 -d'='`
        
        # Get step
        if [ "`echo $line | grep 'span'`" ]
        then
            step=`echo $line | grep -o 'span.*' | cut -f2 -d'='`
        else
            step=1
        fi
    else
        # Writing data into bedgraph file
        if (( $step > 1 ))
        then
            IFS=' ' read position count <<< $line
            for (( i=0; i<$step; i++ ))
            do
                pos=$(($position+$i))
                
                # Mappability
                #echo "$pos $count" | mawk '{OFS="\t"; print $1-1,$1,$2/100}' | sed s/^/"$chromName\t"/ >> $outputFile
                
                # Normal wig file
                echo "$pos $count" | mawk '{OFS="\t"; print $1-1,$1,$2}' | sed s/^/"$chromName\t"/  >> $outputFile
            done   
        else
            # Mappability
            #echo $line | mawk '{OFS="\t"; print $1-1,$1,$2/100}' | sed s/^/"$chromName\t"/ >> $outputFile
            
            # Normal wig file
            echo $line | mawk '{OFS="\t"; print $1-1,$1,$2}' | sed s/^/"$chromName\t"/ >> $outputFile
        fi
    fi
done < $wigFile
#done

EndTime=$(date +%s)
ElapsedTime=$(($EndTime - $StartTime))

echo "File" `basename $wigFile` "converted into $outputFile in $(($ElapsedTime / 60)) minutes and $(($ElapsedTime % 60)) seconds !"

