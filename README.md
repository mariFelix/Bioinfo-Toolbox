# BIOINFO TOOLBOX

2016-11-16

Marianne S. Felix

marianne.sabourin-felix.1@ulaval.ca
------------------------------------

# INTRODUCTION

Various scripts.

# WARNINGS

These scripts are provided without warranty. You can redistribute and/or modify it under the terms of the GNU General Public License v3.


# binningBedGraph.sh

## INTRODUCTION

This script will do the binning (cut the coverage into windows and compute the average for each window) of bedGraph files.

## HOW TO USE IT

```
./binningBedGraph.sh input window output
```
###### EXAMPLE

```
./binningBedGraph.sh UBF_Rep1-MmrDNA+35.bedGraph 25 UBF_Rep1-MmrDNA+35_bin25.bedGraph
```

## HOW IT WORKS



# dataMining_v2.py

## INTRODUCTION

This script allows to screen the GEO datasets database and the SRA database.

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


# sortPeakByFE_v2.py

## INTRODUCTION

## HOW TO USE IT

## HOW IT WORKS


# wigToBedgraph_v2.0.py

## INTRODUCTION

## HOW TO USE IT

## HOW IT WORKS

