# BIOINFO TOOLBOX

2016-11-16

Marianne S. Felix

marianne.sabourin-felix.1@ulaval.ca
------------------------------------

# INTRODUCTION

Various scripts.

# WARNINGS

These scripts are provided without warranty. You can redistribute and/or modify it under the terms of the GNU General Public License v3.

# bamSplit.sh

## INTRODUCTION

This scripts split BAM files by chromosomes.

## DEPENDECIES

* bamtools

## HOW TO USE IT

```
./bamSplit.sh file.bam outputDirectory
```

## HOW IT WORKS

This script uses bamtools to split the input file into chromsomes. It will create files similar to : file-REF_*.bam

# bamToBigWig.sh

## INTRODUCTION

This script converts BAM files into bigWig.

## DEPENDECIES

* BedTools
* bedGraphToBigWig from UCSC (http://hgdownload.soe.ucsc.edu/admin/exe/)

## HOW TO USE IT

```
./bamToBigWig.sh file path/to/bedgraphToBigWigFolder
```

###### EXAMPLE

```
./bamToBigWig.sh file.bam /App/UCSC
```

## HOW IT WORKS

This script first converts the input BAM file into bedGraph and then converts the resulting bedGraph file into bigWig.


# binningBedGraph.sh

## INTRODUCTION

This script will do the binning (cut the coverage into windows and compute the average for each window) of bedGraph files.

## DEPENDECIES

* bc
* mathlib

## HOW TO USE IT

```
./binningBedGraph.sh input window output
```
###### EXAMPLE

```
./binningBedGraph.sh UBF_Rep1-MmrDNA+35.bedGraph 25 UBF_Rep1-MmrDNA+35_bin25.bedGraph
```

## HOW IT WORKS

This script will read the input file line by line and extract the following variables : chromosome, start, end, coverage. It will then compute the mean in each window and replace the coverage by the mean computed. The result will be written in the output file.

# dataMining_v2.py

## INTRODUCTION

This script allows to screen the GEO datasets database and the SRA database.

## DEPENDECIES

* Python
* BioPython

## HOW TO USE IT


```
python dataMining.py
```

```
===============
    M E N U 
  Data Mining 
===============

Please make a choice among the following :
    1 - Search datasets on database (ChIP-seq by protein)
    2 - Search datasets on database (ChIP-seq by cell type)
    3 - Search datasets on database (DNase-seq by cell type)
    4 - Download datasets from tsv file
    5 - Quit the program

Your choice : 
```

###### OPTION 1

```
 1 - Search datasets on database (ChIP-seq by protein)
----------------------------------------------------------
Please enter the maximum number of Id you want to retrieve in each of the 2 databases [1-100000] :
```

```
Please enter a protein name (or many coma-separated) you want to search datasets for :
```

```
Number of datasets found in GEO database : X
Number of datasets found in SRA database : Y
There is Z non overlapping datasets found for this query !
Do you want to save the information of this search in a tsv file ? [yes|no] : 
```

```
Please enter an output file name for the Z dataset(s) found :
```

```
Information extracted for GEO Id : XXXXXXXXX
...
Information extracted for SRA Id : XXXXXX
...
File test.tsv created !
```

###### OPTION 2

```
 2 - Search datasets on database (ChIP-seq by cell type)
----------------------------------------------------------
Please enter the maximum number of Id you want to retrieve in each of the 2 databases [1-100000] :
```

```
Please enter a cell type (or many coma-separated) you want to search datasets for :
```

```
Number of datasets found in GEO database : X
Number of datasets found in SRA database : Y
There is Z non overlapping datasets found for this query !
Do you want to save the information of this search in a tsv file ? [yes|no] :
```

```
Please enter an output file name for the Z dataset(s) found :
```

```
Information extracted for GEO Id : XXXXXXXXX
...
Information extracted for SRA Id : XXXXXXX
...
File test.tsv created !
```

###### OPTION 3

```
 3 - Search datasets on database (DNase-seq by cell type)
----------------------------------------------------------
Please enter the maximum number of Id you want to retrieve in each of the 2 databases [1-100000] :
```

```
Please enter a cell type (or many coma-separated) you want to search datasets for :
```

```
Number of datasets found in GEO database : X
Number of datasets found in SRA database : Y
There is Z non overlapping datasets found for this query !
Do you want to save the information of this search in a tsv file ? [yes|no] :
```

```
Please enter an output file name for the Z dataset(s) found :
```

```
Information extracted for GEO Id : XXXXXXXXX
...
Information extracted for SRA Id : XXXXXXX
...
File test.tsv created !
```

###### OPTION 4

```
 4 - Download datasets from tsv file
----------------------------------------------------------
Download datasets from tsv file
Please enter a tsv file you want to download datasets for :
```

```
Processing GEO Id # XXXXXXXXX
Read Y spots total
Written Y spots total
File XXXXXXXXX.sra downloaded !
...
Processing SRA Id # XXXXXXX
Read Y spots total
Written Y spots total
File XXXXXXXXX.sra downloaded !
...
```

###### OPTION 5

```
Exiting program... Goodbye ! :)
```

## HOW IT WORKS

This script use queries (protein name, cell type, experiment type) to screen into the GEO datasets database and the SRA database of ncbi (https://www.ncbi.nlm.nih.gov/). If a dataset is common to de two databases, only the GEO Id of this dataset will be kept. The resulting datasets found can be written in an output tsv file.

The output contains :
* The database where the data were found
* The Id
* The descriptor
* The format of the data found in the database

If the 4th option is chosen, the user can download the data previously found with option 1 to 3. For this, the script must be rerun with option 4 after being run with option 1 to 3.


# gcContent_Ekblom.sh

## INTRODUCTION

This script is used to get the gc content pattern (% GC) in bedgraph of a sequence using an averaging window set by the user. The input file must contain the postion along the sequence and the corresponding base (A, T, G or C), one position per line. The chromosome name is hard coded but can be changed.

## HOW TO USE IT

```
./gcContent_Ekblom.sh input window output
```

The input file must be in this format :

Position    Base

###### EXAMPLE

```
1   A
2   C
3   C
4   G
5   T
```

In order to have this format, use the following command :

```
sed '1d' rDNA.fa | tr -d '\n' | fold -sw 1 | nl -n ln > rDNA_twoColumns.fa
```

## HOW IT WORKS

This script assign a score according to the base found at every position using a non-overlapping sliding window. G, C and S (G or C) is 1 (100% GC). B (C or G or T) and V (A or C or G) is 0.66666. K (G or T), M (A or C), R (A or G) and Y (C or T) is 0.5. D (A or G or T) and H (A or C or T) is 0.33333. Finally, A, T, W (A or T) and N (A or C or G or T) is 0. The script sums up the scores until reaching the window size. Then it divides the total score of the window by the size of the window and multiply it by 100. It then print the chromosome, the start and the end position of the window and the GC content in a BedGraph format. You can then visualize the resulting file in IGV.


# GCcontentOfPeaksCalledByMACS2.sh

## INTRODUCTION

This script computes the average GC% for each peak regions called by MACS2.

## HOW TO USE IT

This script requires the genome files (*.fa) to be in one line. To do so, use the following command :

```
sed '1d' file.fa | tr -d '\n' > output.fa
```

###### OR

```
for file in *.fa; do sed '1d' $file | tr -d '\n' > outputDirectory/$file; done
```


## HOW IT WORKS

This script reads the *.narrowPeak file line by line (one peak per line) and extracts the sequence for each peak from the reference sequence file. It then computes the average GC% for this sequence. The result is printed in the output file along with the chromosome name, the start and the end positions of the peak, the peak name, its sequence and the computed GC%. At the end, the average GC% of all the peaks is printed to the terminal.


# normalizationWithqPCR.py

## INTRODUCTION

This script allows to normalize bedgraphs KO files according to the qPCR results of mouse rDNA primers (IGS3, SpPr, Tsp, Prom, 47S, ETS, 28S, T1A).

## HOW TO USE IT

```
python normalizationWithqPCR.py
```

###### Enter the input files

```
Please enter your non-induced normalized file :
```

```
Please enter your induced normalized file :
```

```
Please enter your output file name :
```

###### Enter the values of qPCR for the differents primers you want to use

```
Do you want to use the IGS3 primer [yes|no] :
```

```
Please enter the qPCR value for IGS3 non-induced :
```

```
Please enter the qPCR value for IGS3 induced :
```

```
...
```

```
Do you want to use the TIA primer [yes|no] :
```

```
Please enter the qPCR value for TIA non-induced :
```

```
Please enter the qPCR value for TIA induced :
```

###### The final ratio is printed

```
Final mean : X
```

## HOW IT WORKS

This script ask the user for the qPCR values of non-induced and induced (KO) datasets for different primers on the mouse rDNA. One can choose any primer(s) to do the analysis. The script stocks the ratio between the non-induced over the induced in a dictionnay of lists with 2 argument along with the position of the primers (start of the forward primer to the start of the reverse primer). The primers with a ratio of zero (not chosen by the user) are removed from the dictionnary. Then the script reads the non-induced input file and the induced input file line by line and for each primer, it stocks the coverage that falls inside the primer region do to an averaging. The average coverage of the non-induced is then divided by the average coverage of the induced to create a ratio. The ratio of the ChIP-seq datasets then created is divided by the ratio of the qPCR computed with the input values for each primer. The "final mean" is then generated by doing an averaging of all the ratio ChIP/qPCR. Finally, the induced file is multiplied by the final mean and printed in the output file spectified by the user.

# removePrimerContamination.sh

## INTRODUCTION

This script removes the contamination of primers on the rDNA in the aligned BAM files.

## HOW TO USE IT

```
./removePrimerContamination.sh processedFile1.bam processedFile2.bam processedFile3.bam processedFile4.bam outputFolder
```

## HOW IT WORKS

This script takes as an input 4 files (hard coded) and removes the forward and reverse primer sequences (sequences are hard coded but can be modified according to the primer contamining the datasets). The BAM datasets are converted to SAM format. Then a "grep -v" is done to extract reads that DOESN'T contain the forward and reverse primer sequences. The result is then printed in the ouput file (in SAM format) that is reconverted into BAM format.

# sortPeakByFE_v2.py

## INTRODUCTION

This script allows to sort peaks according to their fold enrichment (FE).

## HOW TO USE IT

```
$ python sortPeakByFE.py --help

usage: sortPeakByFE.py [-h] (-r | -k) threshold input output

  This script allows to sort peaks according to their fold enrichment (FE).

positional arguments:
  threshold   Fold enrichment threshold (float over 1)
  input       Input file
  output      Output file (will be overwritten if exists)

optional arguments:
  -h, --help  show this help message and exit
  -r          Removes peaks with FE below threshold
  -k          Keeps peaks with FE of IP over (threshold * FE of KO)

examples:
  python sortNarrowPeakByFE.py -r threshold in.narrowPeak out.narrowPeak
       * "in.narrowPeak" is the standard output from MACS2 callpeak

  python sortNarrowPeakByFE.py -k threshold in.narrowPeak out.narrowPeak
       * "in.narrowPeak" is the output from bedtools intersect -loj
```

###### -r option (removeFPnarrowPeak)

This function removes peaks that have a fold enrichment below a given threshold. It prints the resulting peaks into the outputFile and returns the number of peak over the given threshold.

###### -k option (keepFNnarrowPeak)

This function keeps peaks that have a fold enrichment over (threshold * fold enrichment of KO). It prints the resulting peaks into the outputFile and returns the number of peak over (threshold * fold enrichment of KO).

## HOW IT WORKS

###### With the -r option

This script takes as an input a threshold, a standard narrowPeak file and an ouput file name. It stores the lines of the input file and sort them by fold enrichment (FE, 6th column). It then kepps only the peaks with an FE over the given threshold and output them into the ouput file.

###### With the -k option

*** Under construction ***

# wigToBedgraph_v2.0.py

## INTRODUCTION

This script converts *.wig files into *.bedgraph files.

## HOW TO USE IT

```
./wigToBedgraph.sh inFile.wig outFileName
```

###### EXAMPLE

```
./wigToBedgraph.sh MmrDNA35.wig experiment1
```

## HOW IT WORKS

This script reads the *.wig file line by line. It extracts the chromosome name and the span (number of base having the same coverage) and the coverage. Then it prints the position (base by base) along with its coverage in the output *.bedgraph file.
