#!/usr/bin/python
# encoding: utf-8

# normalizationWithqPCR
# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 1.0
# 2018-02-15
#
# [CAUTION] : This script works only if the ratio of the qPCR is greater than the ratio of the ChIP-seq.


"""
normalizationWithqPCR module
"""

import os
import sys
import subprocess


##################################
#                                #
#####        Primers         #####
#                                #
##################################

primerIGS3 = [12158,12414]
primerSpPr = [12581,12784]
primerTsp = [12758,12952]
primerProm = [14639,14854]
primer47S = [14973,15134]
primerETS = [17892,18035]
primer28S = [25029,25225]
primerT1A = [28230,28420]


##################################
#                                #
#####        Functions       ##### 
#                                #
##################################

def primerValuesInput(primerName):

	while True:

		rep = raw_input("Do you want to use the {} primer [yes|no] : ".format(primerName))

		if rep == "yes":
			primerN = raw_input("Please enter the qPCR value for {} non-induced : ".format(primerName))
			primerI = raw_input("Please enter the qPCR value for {} induced : ".format(primerName))
			break
		elif rep == "no":
			break
		print "[ERROR] : The answer must be yes or no !"

	if rep == "yes":
		return float(primerN)/float(primerI)
	else:
		return None



##################################
#                                #
#####     User's input       #####
#                                #
##################################

# Caution
print "[CAUTION] : This script works only if the ratio of the qPCR is greater than the ratio of the ChIP-seq."

# Input and output files
nonInduced = raw_input("Please enter your non-induced normalized file : ")
induced = raw_input("Please enter your induced normalized file : ")
output = raw_input("Please enter your output file name : ")

# Values of IP non-induced and induced at primer location
valIGS3 = primerValuesInput("IGS3")
valSpPr = primerValuesInput("SpPr")
valTsp = primerValuesInput("Tsp")
valProm = primerValuesInput("Prom")
val47S = primerValuesInput("47S")
valETS = primerValuesInput("ETS")
val28S = primerValuesInput("28S")
valT1A = primerValuesInput("TIA")

# Creation of a dictionnary
dictOfPrimers = {"IGS3" : [primerIGS3, valIGS3], "SpPr" : [primerSpPr, valSpPr],
				 "Tsp" : [primerTsp, valTsp], "Prom" : [primerProm, valProm], 
				 "47S" : [primer47S, val47S], "ETS" : [primerETS, valETS], 
				 "28S" : [primer28S, val28S], "T1A" : [primerT1A, valT1A]}

# Remove empty qPCR value
for x in dictOfPrimers.keys():
	if dictOfPrimers[x][1] == None:
		del dictOfPrimers[x]


##################################
#                                #
#####        M A I N         #####
#                                #
##################################

# Variable initiation
meanTotalN = []
meanTotalI = []
dictChIPratio = {}

# Open non-induced and induced files
with open(nonInduced ,"r") as n, open(induced, "r") as i:

	# Stock lines into variables
	linesN = n.readlines()
	linesI = i.readlines()

	# Iterate through the dictionnary
	for x in dictOfPrimers.keys():

		start = dictOfPrimers[x][0][0]
		end = dictOfPrimers[x][0][1]

		# Variable initiation, reinitiate after each primer
		arrN = []
		arrI = []

		# Stock the coverage inside the primer position to do the mean
		for y in range(int(start), int(end+1), 1):

			# Assign values
			chromN, posStartN, posEndN, covN = linesN[y].split()
			chromI, posStartI, posEndI, covI = linesI[y].split()

			# Append the coverage to a list
			arrN.append(covN)
			arrN = [float(i) for i in arrN]

			arrI.append(covI)
			arrI = [float(i) for i in arrI]

		# Compute the mean of the coverage for a primer
		meanN = sum(arrN) / float(len(arrN))
		meanI = sum(arrI) / float(len(arrI))

		# Compute the ratio of non-induced/induced
		ratioChIP = meanN/meanI
		
		# Append the primer and the ratio to a list
		dictChIPratio[x] = ratioChIP 

	# Variable initiation
	listRatio = []

	for z in dictChIPratio.keys():

		# Compute the ration of ChIP-seq/qPCR at each primer
		chip = dictChIPratio[z]
		qPCR = dictOfPrimers[z][1]
		ratio = float(chip) / float(qPCR)

		#print "Mean ChIP-seq for {} : {}".format(z, chip)
		#print "Mean qPCR for {} : {}".format(z, qPCR)
		#print "Ratio for {} : {}".format(z, ratio)

		# Append the final ratio to a list
		listRatio.append(ratio)

	# Compute the mean of all the ratio ChIP-seq/qPCR for all the primers
	finalMean = sum(listRatio) / float(len(listRatio))

	print "Final mean :", finalMean

# Writing output file (coverage * finalMean)
with open(induced, "r") as i, open(output, "a") as o:

	# Extract header
	header = i.readline()

	# Write header to output file
	o.write(header)

	# Read through the induced file line by line
	for line in i:

		# Assign values
		chrom, start, end, cov  = line.split()

		# Compute new coverage
		normCov = float(cov) * float(finalMean)

		# Write line to output file
		o.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, normCov))









