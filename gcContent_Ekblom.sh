#!/bin/bash

# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 1.0
# 2016-05-05
#
# gcContent_Ekblom.sh

input=$1
window=$2
output=$3
chr="MmrDNA+35"

echo -e "track type=bedgraph visibility=full color=0,0,204 altColor=204,0,0 autoScale=on maxHeightPixels=128:128:11 viewLimits=0.0:25.0 yLineMark=0.0 windowingFunction=maximum" > $output

i=0
while read -r line;
do    
    # Get line arguments
    read pos base <<< $line
    
    #echo "pos : $pos, base : $base"
    
    case $base in
    
        C|c|G|g|S|s)
            # S -> C or G
            count=1
            ;;
        
        B|b|V|v)
            # B -> C or G or T
            # V -> A or C or G
            count=0.66666
            ;;
        
        K|k|M|m|R|r|Y|y)
            # K -> G or T
            # M -> A or C
            # R -> A or G
            # Y -> C or T
            # N -> A or C or G or T
            count=0.5
            ;;
        
        D|d|H|h)
            # D -> A or G or T
            # H -> A or C or T
            count=0.33333
            ;;
        
        A|a|T|t|W|w|N|n)
            # W -> A or T
            # N -> A or C or G or T
            count=0
            ;; 
        
        *)
            echo "`tput setaf 1`[ERROR101]`tput sgr0` : An internal error occured, please report it to marianne.sabourin-felix.1@ulaval.ca"
            exit
            ;;
    esac
    
    sumCount=$((sumCount+count))
    
    i=$((i+1))
    
    if [[ $i == $window ]]
    then
        gcContent=`bc -l <<< "$sumCount/$window*100"`
        
        echo "$chr $pos $window $gcContent" | mawk '{OFS="\t"; print($1,$2-$3,$2,$4)}' >> $output
        
        i=0
        sumCount=0
        gcContent=0
    fi
done < $input



