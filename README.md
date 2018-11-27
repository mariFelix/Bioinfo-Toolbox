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

If the 4th option is chosen, the user can download the data previously found with option 1 to 3. For this, the script must be rerun with option 4 after being run wiht option 1 to 3.

*test

# extendSmooth.py

## INTRODUCTION

## HOW TO USE IT

## HOW IT WORKS


# gcContent_Ekblom.sh

## INTRODUCTION

## HOW TO USE IT

## HOW IT WORKS


# normalizationWithqPCR.py

## INTRODUCTION

## HOW TO USE IT

## HOW IT WORKS


# removePrimerContamination.sh

## INTRODUCTION

## HOW TO USE IT

## HOW IT WORKS


# sortPeakByFE_v2.py

## INTRODUCTION

## HOW TO USE IT

## HOW IT WORKS


# wigToBedgraph_v2.0.py

## INTRODUCTION

## HOW TO USE IT

## HOW IT WORKS

