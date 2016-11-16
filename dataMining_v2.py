#!/usr/bin/python
# encoding: utf-8

# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 2.0
# 2016-05-30
#
# dataMining.py
# Tested on Ubuntu 14.04 LTS with Python 2.7.11
################################################################################

"""
Data mining module
"""

import sys
from subprocess import call
from Bio import Entrez
#import time
import re
from urllib2 import urlopen


class color:
    """
    Print in color
    """
    bold = '\033[1m'
    red = '\033[91m'
    end = '\033[0m'


def exceptionHandling(valueType, value, traceback):
    """
    Intern function that catches unhandled exceptions.
    """
    errorMess = ("[ERROR] An internal error occured, please contact "
                 "marianne.sabourin-felix.1@ulaval.ca and mention the "
                 "following DEBUG message. Thank you !")

    # Formatted string showing the error type
    formType = str(valueType).strip('>').split(' ')[1].strip("'")

    print errorMess
    print "[DEBUG] {} : {}".format(formType, value)


def _getOrgnByCellType(cellType):
    """
    Intern function that returns the organism
    corresponding to the cell type provided.

    Return organism abbreviation
    """

    # List of cell type by organism
    #TODO : Complete lists
    #TODO : Case insensitive (k562 == K562)
    mm = ['MEF', 'MEL', 'CH12', 'NIH3T3', 'mES']
    hs = ['K562', 'HEK293T', 'HMEC', 'HMLER']

    if cellType in mm:
        Orgn = "Mus+musculus[Orgn]"
    elif cellType in hs:
        Orgn = "Homo+sapiens[Orgn]"
    else:
        Orgn = ""

    return Orgn


def _getSeqMeth(seqMeth):
    """
    Intern function that returns
    the sequencing method

    Return the sequencing method
    """

    if seqMeth == "chipSeq":
        seq = "ChIP-seq OR ChIPseq OR chip+seq OR ChIP+sequencing"
    elif seqMeth == "dnaseSeq":
        seq = "DNAse1 OR DNAse-seq OR DNAseSeq OR dnase+seq OR DNAse+sequencing"

    return seq


def _datasetInfoGEO(Id):
    """
    This function returns a formated line for a given GEO Id

    Return Database Id Title DataFormat
    """

    summaryTmp = Entrez.esummary(db="gds", Id="{}".format(Id))
    summary = Entrez.read(summaryTmp)

    # UnicodeEncodeError handled by *.encode('ascii', 'replace')
    return "GeoDatasets\t{0}\t{1}\t{2}\n".format(Id,
            summary[0]["title"].encode('ascii', 'replace'),
            summary[0]["suppFile"].encode('ascii', 'replace'))


def _datasetInfoSRA(Id):
    """
    This function returns a formated line for a given SRA Id

    Return Database Id Title
    """

    summaryTmp = Entrez.esummary(db="sra", Id="{}".format(Id))
    summary = Entrez.read(summaryTmp)

    parsing = summary[0]["ExpXml"]
    title = parsing[parsing.find("<Title>") + 7:parsing.find("</Title>")]
    #formattedTitle = title.split(':', 1)[1].encode('ascii', 'replace')

    # UnicodeEncodeError handled by *.encode('ascii', 'replace')
    return "SRA\t{0}\t{1}\tSRA\n".format(Id, title.encode('ascii', 'replace'))


def getUserInput():
    """
    Function that prints the main menu

    Return user choice
    """

    mainMenu = (color.bold + "=" * 15 + "\n    M E N U \n  Data Mining \n"
                + "=" * 15 + color.end +
                "\n\nPlease make a choice among the following :")

    menuChoices = {
        "1": " 1 - Search datasets on database (ChIP-seq by protein)",
        "2": " 2 - Search datasets on database (ChIP-seq by cell type)",
        "3": " 3 - Search datasets on database (DNase-seq by cell type)",
        "4": " 4 - Download datasets from tsv file",
        "5": " 5 - Quit the program"}

    print mainMenu

    # Print choices from dictonary
    for choice in sorted(menuChoices):
        print " " * 3 + menuChoices[choice]

    # Validation of user input
    choice = raw_input("\nYour choice : ")
    while choice not in menuChoices:
        print(color.red + "[ERROR] " + color.end +
              "Your choice must be an integer between 1 and 5 !")
        choice = raw_input("Your choice : ")

    # Print header of choice
    if choice == "5":
        print "\nExiting program... Goodbye ! :)"
        sys.exit()
    else:
        print("\n" + color.bold + menuChoices[choice] + "\n"
              + "-" * 58 + color.end)

    return choice


def searchChIPSampleGEObyProt(protein, retMax):
    """
    This function search for ChIP-seq Geo Sample in the GEO database for
    Mus musculus and Homo sapiens for a given protein.

    Return GeoSamples Id
    """

    #Entrez.email = "marianne.sabourin-felix.1@ulaval.ca"
    seq = "(ChIP-seq OR ChIPseq OR chip+seq OR ChIP+sequencing)"
    orgn = "Mus+musculus[Orgn] OR Homo+sapiens[Orgn]"

    handle = Entrez.esearch(db="gds", retmode="xml", retmax=retMax,
             term="gsm[Filter] AND ({}) AND ({}) AND ({})"
             .format(seq, orgn, protein))
    record = Entrez.read(handle)

    return record["IdList"]


def searchSampleGEObyCellType(cellType, seqMeth, retMax):
    """
    This function search for Geo Sample in the GEO database for
    Mus musculus and/or Homo sapiens for a given cell type.

    Return GeoSamples Id
    """

    # Get sequencing method
    seq = _getSeqMeth(seqMeth)

    # Format cellType
    splitCell = cellType.split(' OR ')
    sourceName = [cell + "[Source Name]" for cell in splitCell]
    formattedCell = ' OR '.join(sourceName)

    # Get organism name(s)
    # Portable for many organism but maybe too complicated for the actual code
    #TODO : If one cell return "" and the other not, put mm and hs...
    notFound = True
    setOfCellTypes = set()
    for cell in splitCell:
        knownOrgn = _getOrgnByCellType(cell)
        if len(knownOrgn) != 0:
            setOfCellTypes.add(knownOrgn)
            notFound = False

    # Format organism name(s)
    if notFound is True:
        orgn = "Mus+musculus[Orgn] OR Homo+sapiens[Orgn]"
    else:
        orgn = ' OR '.join(setOfCellTypes)

    # Search into database
    handle = Entrez.esearch(db="gds", retmode="xml", retmax=retMax, term=
    "gsm[Filter] AND ({}) AND ({}) AND ({})".format(seq, orgn, formattedCell))

    record = Entrez.read(handle)

    return record["IdList"]


def searchSampleSRA(query, seqMeth, retMax):
    """
    This function searches for ChIP-seq or DNase-seq sample in the SRA database
    for Mus musculus and/or Homo sapiens for a given protein or cell type.

    Return a list of SRA Id
    """

    # Get sequencing method
    seq = _getSeqMeth(seqMeth)

    #TODO : Put in a function
    # If query is a cell type, get organism name
    # Portable for many organism but maybe too complicated for the actual code
    splitCell = query.split(' OR ')
    notFound = True
    setOfCellTypes = set()
    for cell in splitCell:
        knownOrgn = _getOrgnByCellType(cell)
        if len(knownOrgn) != 0:
            setOfCellTypes.add(knownOrgn)
            notFound = False

    # Format organism name(s)
    if notFound is True:
        orgn = "Mus+musculus[Orgn] OR Homo+sapiens[Orgn]"
    else:
        orgn = ' OR '.join(setOfCellTypes)

    # Fetch list of SRA Id
    handle = Entrez.esearch(db="sra", retmode="xml", retmax=retMax,
             term="({}) AND ({}) AND ({})".format(seq, orgn, query))
    record = Entrez.read(handle)

    return record["IdList"]


def clearIdList(listGEO, listSRA):
    """
    This function removes from the SRA list Id that is also in the GEO list.
    (Avoid doing the whole analysis if already present in GEO)

    Return a list of SRA Id (minus the ones with GEO datasets)
    """

    # Stock GSM ID
    listOfGSM = []
    for geoDataset in listGEO:

        summaryTmp = Entrez.esummary(db="gds", Id="{}".format(geoDataset))
        summary = Entrez.read(summaryTmp)

        listOfGSM.append(summary[0]["Accession"])

    # Clean SRA list
    alreadyInGEO = []
    for sraDataset in listSRA:

        summaryTmp = Entrez.esummary(db="sra", Id="{}".format(sraDataset))
        summary = Entrez.read(summaryTmp)

        accession = summary[0]["ExpXml"]

        for GSM in listOfGSM:
            if GSM in accession:
                alreadyInGEO.append(sraDataset)

    return [x for x in listSRA if x not in alreadyInGEO]


def createTsvList(listGEOid, listSRAid, outputFile):
    """
    This function creates a tsv file of GEO and SRA Id
    for a particular search.

    Return ???
    """

    with open(outputFile, 'w') as o:
        for GEOid in listGEOid:
            formattedInfoGEO = _datasetInfoGEO(GEOid)
            o.write(formattedInfoGEO)
            print "Information extracted for GEO Id : {}".format(GEOid)

        for SRAid in listSRAid:
            formattedInfoSRA = _datasetInfoSRA(SRAid)
            o.write(formattedInfoSRA)
            print "Information extracted for SRA Id : {}".format(SRAid)

    print "File {} created !".format(outputFile)


def fetchGEOdata(Id):
    """
    This function fetch data from a given GEOid

    Download data from FTP server
    """

    # Fetch ftpLink from GEOid
    summaryTmp = Entrez.esummary(db="gds", Id="{}".format(Id))
    summary = Entrez.read(summaryTmp)
    ftpLink = summary[0]["FTPLink"]

    # Variable initialisation
    listOfFiles = []
    newFtpLink = []
    isAfolder = True

    # Fetch all files in every directories
    while isAfolder is True:
        isAfolder = False

        # Convert ftpLink to list if not already
        if type(ftpLink) is not list:
            ftpLink = ftpLink.split()

        for link in ftpLink:
            # Look into the ftpLink to list the files and/or folders
            urlpath = urlopen(link)
            ls = urlpath.read()
            folderLines = ls.splitlines()

            for lines in folderLines:
                # Fetch the file or folder name
                elemName = lines.split()[-1]

                # If elemName is a folder, add to newFtpLink
                if lines[0] == 'd':
                    newFtpLink.append(link + elemName + "/")
                    isAfolder = True
                # If elemName is a file, add to listOfFiles
                else:
                    listOfFiles.append((link, elemName))

        # Assign the new ftpLink list and clear the variable
        if isAfolder is True:
            ftpLink = newFtpLink[:]
            del newFtpLink[:]

    # Download files
    for tup in listOfFiles:
        ftpLink = tup[0] + tup[1]
        output = urlopen(ftpLink)

        with open(tup[1], 'wb') as reading:
            reading.write(output.read())

        print "File {} downloaded !".format(tup[1])

    # Find a other return value (failure of success ?!)
    return listOfFiles


def returnSRAaccession(Id):
    """
    This function return SRA accession number(s) for a given Id

    Return a list of SRA accession number(s)
    """

    #TODO: Find a more elegant way to do it...

    summaryTmp = Entrez.esummary(db="sra", Id="{}".format(Id))
    summary = Entrez.read(summaryTmp)

    parsing = summary[0]["Runs"]

    listOfIndex = [index.start() for index in re.finditer("<Run acc=", parsing)]

    listOfSRAaccession = []
    for index in listOfIndex:
        sra = parsing[index:len(parsing)].split("\"")[1]
        listOfSRAaccession.append(sra)

    return listOfSRAaccession


def main():
    """
    Main function
    """
    #TODO : Create a function to get user email ?
    Entrez.email = "marianne.sabourin-felix.1@ulaval.ca"

    #Id = raw_input("SRA Id : ")
    #listSRA = returnSRAaccession(Id)
    #for sra in listSRA:
        #print "SRA for Id {} : {}".format(Id, sra)

    #sys.exit()

    # Catch unhandled exceptions ### ADD IF RELEASED
    #sys.excepthook = exceptionHandling

    #TODO : Create a function to get user email ?
    Entrez.email = "marianne.sabourin-felix.1@ulaval.ca"

    # Get user input
    choice = getUserInput()

    # Search datasets on database
    if choice != "4":
        # Get sequencing method
        seq = ("dnaseSeq" if choice == "3" else "chipSeq")

        #TODO : Add validation of user's input
        # Get the maximum number of Id to retrieve (retmax)
        retMax = raw_input("Please enter the maximum number of Id you want to "
                           "retrieve in each of the 2 databases [1-100000] : ")

        # Search by protein name
        if choice == "1":
            # Get protein(s) name(s)
            query = raw_input("Please enter a protein name "
                              "(or many coma-separated) you want to search "
                              "datasets for : ").split(',')

        # Search by cell type
        elif choice in ("2", "3"):
            # Get cell type(s)
            query = raw_input("Please enter a cell type "
                              "(or many coma-separated) you want to search "
                              "datasets for : ").split(',')

        # Format user input
        if len(query) == 1:
            queryFormatted = ' '.join(query)
        else:
            querySpace = [space.strip(' ').replace(' ', '+') for space in query]
            queryFormatted = ' OR '.join(querySpace)

        # Get list GEO Id
        #TODO : Merge the two GEO function into one
        if choice == "1":
            listGEO = searchChIPSampleGEObyProt(queryFormatted, retMax)
        else:
            listGEO = searchSampleGEObyCellType(queryFormatted, seq, retMax)

        # Get list SRA Id
        listSRA = searchSampleSRA(queryFormatted, seq, retMax)

        # Remove from SRA list the data already in GEO list
        purgedListSRA = clearIdList(listGEO, listSRA)

        # Check if the research returns results
        numOfHit = len(listGEO + purgedListSRA)

        if numOfHit == 0:
            print "No dataset found for this query !"

        else:
            # Print the number of hit
            print "Number of datasets found in GEO database : {}" \
                  .format(len(listGEO))
            print "Number of datasets found in SRA database : {}" \
                  .format(len(listSRA))

            print "There is {} non overlapping datasets found for this query !"\
                  .format(numOfHit)

            #for dataset in listGEO + listSRA:
                #print dataset

            # Save information or not
            save = raw_input("Do you want to save the information of this "
                             "search in a tsv file ? [yes|no] : ")

            while save not in ("yes", "no"):
                print(color.red + "[ERROR] " + color.end +
                      "Your choice must be yes or no !")
                save = raw_input("Do you want to save the information of this "
                                 "search in a tsv file ? [yes|no] : ")

            if save == "yes":
                # Get output filename from user input
                output = raw_input("Please enter an output file name for the "
                                   "{} dataset(s) found : ".format(numOfHit))

                # Create tsv file
                createTsvList(listGEO, purgedListSRA, output)

            else:
                print "The information was not saved !"

    # Download datasets from tsv file
    else:
        #TODO
        print "Download datasets from tsv file"

        # Get input file (tsv file from step A, manually sorted)
        inputFile = raw_input("Please enter a tsv file you want to "
                              "download datasets for : ")

        with open(inputFile, 'r') as inFile:
            for data in inFile:

                # Assign elements
                line = data.split()
                database = line[0]
                Id = line[1]

                # Download datasets
                if database == "GeoDatasets":
                    print "Processing GEO Id # {}".format(Id)
                    fetchGEOdata(Id)

                elif database == "SRA":
                    print "Processing SRA Id # {}".format(Id)
                    listSRAacc = returnSRAaccession(Id)
                    for acc in listSRAacc:
                        call(["fastq-dump", "--gzip", acc, "2> /dev/null"])
                        print "File {}.sra downloaded !".format(acc)

                else:
                    print " [ERROR] : Database not supported !"


if __name__ == '__main__':

    main()
