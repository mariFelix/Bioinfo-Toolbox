#!/bin/bash

# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 1.0
# 2017-10-16
#
# binning25bp.sh

input=$1
window=$2
output=$3

echo -e "track type=bedgraph visibility=full color=0,0,204 altColor=0,0,204 autoScale=on maxHeightPixels=128:128:11 viewLimits=0.0:25.0 yLineMark=0.0 windowingFunction=maximum" > $output

i=1
count=0
windowNew=$window
while read -r line;
do    
    # Get line arguments
    IFS=' ' read chr start end value <<< $line
    
    #echo "chr : $chr, start : $start, end : $end, value : $value"

    count=`bc -l <<< "$count+$value"`
    
    if [[ $windowNew == $end ]]
    then

        mean=`bc -l <<< "$count/$window"`
        
        count=0
        
        i=$((i+1))
        windowNew=$((window*i))
        
        echo "$chr $end $window $mean" | mawk '{OFS="\t"; print($1,$2-$3,$2,$4)}' >> $output
        
    fi

done < $input
