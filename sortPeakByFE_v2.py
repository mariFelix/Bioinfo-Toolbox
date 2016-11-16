#!/usr/bin/python
# encoding: utf-8

# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# version : 2.0
# 2016-06-10
#
# sortPeakByFE.py
# Tested on Ubuntu 14.04 LTS with Python 2.7.11 adn Python 3.4.3

#TODO : Find a synonym of sort, decrease, extract...

"""
sortPeakByFE module
"""

import os
import sys
import time
import argparse


def printMenu():
    """
    This function print the script menu
    """

    menu = ("\nsortNarrowPeakByFE.py" "\n-------------------------\n"
    "This script allows to sort peaks according to their fold enrichment (FE)."
    "\n\n[USAGE] python sortNarrowPeakByFE.py -option in.narrowPeak "
    "out.narrowPeak threshold\n\n"
    "[OPTIONS] -h/--help : Print this menu\n"
    "                 -r : Removes peaks with FE below threshold\n"
    "                 -k : Keeps peaks with FE of IP over "
                          "(threshold * FE of KO)\n\n"

    "[EXAMPLES] python sortNarrowPeakByFE.py -r in.narrowPeak out.narrowPeak "
    "threshold\n               * \"in.narrowPeak\" is the standard output "
                               "from MACS2 callpeak\n\n"
    "           python sortNarrowPeakByFE.py -k in.narrowPeak out.narrowPeak "
    "threshold\n               * \"in.narrowPeak\" is the output from bedtools"
                               " intersect -loj\n")

    print(menu)
    sys.exit()


def _threshold(value):
    """
    Intern function that raises an error if
    threshold is below or equal to 1

    Return threshold
    """

    try:
        value = float(value)
    except:
        raise argparse.ArgumentTypeError("The threshold must be a float.")

    if value <= 1:
        raise argparse.ArgumentTypeError("The threshold must be over 1.")

    return value


def _validFile(path):
    """
    Intern function that raises an error if filename
    doesn't exists or is not in narrowPeak format

    Return path
    """

    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError("The input file doesn't exist.")
    else:
        try:
            open(path, 'r')
        except:
            raise argparse.ArgumentTypeError("The input file can't be opened"
                                             " in reading mode.")

        if not path.endswith(".narrowPeak"):
            raise argparse.ArgumentTypeError("The input file must be in "
                                             "narrowPeak format.")

    #TODO : Count the number of column in the input file ?

    return path


def parseArgv(argv=None):
    """
    This function allows to parse the command line input options

    Return parsed input parameters
    """

    # Parse arguments
    parser = argparse.ArgumentParser(formatter_class=
                                     argparse.RawDescriptionHelpFormatter,
        description="  This script allows to sort peaks according to their fold"
                  + " enrichment (FE).",
        epilog="""examples:
  python sortNarrowPeakByFE.py -r threshold in.narrowPeak out.narrowPeak
       * "in.narrowPeak" is the standard output from MACS2 callpeak

  python sortNarrowPeakByFE.py -k threshold in.narrowPeak out.narrowPeak
       * "in.narrowPeak" is the output from bedtools intersect -loj""")

    # Mutually exclusives options
    options = parser.add_mutually_exclusive_group(required=True)
    options.add_argument("-r", action="store_true",
                        help="Removes peaks with FE below threshold")
    options.add_argument("-k", action="store_true",
                        help="Keeps peaks with FE of IP over "
                           + "(threshold * FE of KO)")

    parser.add_argument("threshold", type=_threshold,
                        help="Fold enrichment threshold (float over 1)")
    # Positional arguments
    parser.add_argument("input", type=_validFile, help="Input file")
    parser.add_argument("output", type=str, help="Output file (will be "
                                               + "overwritten if exists)")

    return parser.parse_args(argv)


def removeFPnarrowPeak(inputFile, outputFile, threshold):
    """
    This function removes peaks that have a
    fold enrichment below a given threshold

    Print resulting peaks to outputFile
    Return number of peak over the threshold
    """

    # Create a list of list from each line of the input file
    with open(inputFile, 'r') as f:
        lines = [line.split() for line in f]

    # Sort lines by descending fold enrichment
    sortedLines = sorted(lines, key=lambda i: float(i[6]), reverse=True)

    peaks = 0
    with open(outputFile, 'w') as o:
        for line in sortedLines:
            if float(line[6]) < float(threshold):
                break
            else:
                o.write('{}\n'.format('\t'.join(line)))
                peaks += 1

    return peaks


def keepFNnarrowPeak(inputFile, outputFile, threshold):
    """
    This function keeps peaks that have a fold enrichment
    over (threshold * fold enrichment of KO)

    Print resulting peaks to outputFile
    Return number of peak over (threshold * fold enrichment of KO)
    """
    #TODO : Check if file contain more than 10 lines
    #TODO :    (IndexError: list index out of range)

    # Create a list of list from each line of the input file
    with open(inputFile, 'r') as f:
        unique = []
        for line in f:

            # Assignation of FE columns
            elements = line.split()
            FEofIP = float(elements[6])
            # Convert NULL (.) to zero for unique peak (if no KO peak found)
            FEofKO = (float(0) if elements[16] == "." else float(elements[16]))
            multiple = float(threshold) * FEofKO

            # Add left outer join to list if not already in it
            # Peaks in IP may be duplicated (intersect with several KO peaks)
            if FEofIP >= multiple:
                newLine = elements[0:10]
                if newLine not in unique:
                    unique.append(newLine)

    # Sort lines by descending fold enrichment
    sortedUnique = sorted(unique, key=lambda i: float(i[6]), reverse=True)

    # Write to output file
    peaks = 0
    with open(outputFile, 'w') as o:
        for line in sortedUnique:
            o.write('{}\n'.format('\t'.join(line)))
            peaks += 1

    return peaks


def Oldmain():
    """
    Old main function

    """

    # Validation of the number of input parameters
    # TODO : Validation of input parameters
    if len(sys.argv) == 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        printMenu()

    elif len(sys.argv) != 5:
        print((" [ERROR] Incorrect amount of arguments : requires 4, {} \
        provided.".format(len(sys.argv) - 1)))
        sys.exit()

    else:
        # Parameters assignation
        #TODO : find a way for autocompletion or put the option at the end...
        option = sys.argv[1]
        inputFile = sys.argv[2]     # Input file (narrowPeak from MACS2)
        outputFile = sys.argv[3]    # Ouput file (peaks over minFE)
        threshold = sys.argv[4]     # Mininum Fold Enrichment
                                    #(threshold to keep peaks)

    # Function execution according to user's option
    startTime = time.time()

    if option == "-r":
        peaks = removeFPnarrowPeak(inputFile, outputFile, threshold)

    elif option == "-k":
        peaks = keepFNnarrowPeak(inputFile, outputFile, threshold)

    else:
        print((" [ERROR] Option \"{}\" invalid !".format(option)))
        printMenu()

    endTime = time.time()
    elapsedTime = endTime - startTime

    print(("File {} created in {:.5f} seconds with {} peaks found !".
    format(outputFile, elapsedTime, peaks)))


def main():
    """
    Main function

    """

    # Parse argv
    argv = parseArgv()

    # Function execution according to user's option
    startTime = time.time()

    # Call of functions
    if argv.r:
        peaks = removeFPnarrowPeak(argv.input, argv.output, argv.threshold)
    elif argv.k:
        peaks = keepFNnarrowPeak(argv.input, argv.output, argv.threshold)

    endTime = time.time()
    elapsedTime = endTime - startTime

    print(("File {} created in {:.5f} seconds with {} peaks found !".
    format(argv.output, elapsedTime, peaks)))


##################################
#                                #
##########  M A I N  #############
#                                #
##################################

if __name__ == '__main__':

    main()

##################################
