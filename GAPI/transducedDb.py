#!/usr/bin/env python
#coding: utf-8

## DEPENDENCIES ##
# External
import argparse
import os
import sys
import subprocess
import mappy as mp

# Internal
from GAPI import unix
from GAPI import formats


######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('input', help='')
parser.add_argument('index', help='')

parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='Output directory. Default: current working directory')

## 2. Parse user´s input and initialize variables ##
args = parser.parse_args()
inputFile = args.input
index = args.index
outDir = args.outDir


##############################################
## Display configuration to standard output ##
##############################################
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('input: ', inputFile)
print('index: ', index)
print('outDir: ', outDir, "\n")

print('***** Executing ', scriptName, '.... *****', "\n")


##########
## CORE ## 
##########

## 1. Extract sequence for transductions intervals and generate fasta file object
##################################################################################

## Initialize FASTA object
FASTA = formats.FASTA()

## Load Index
index = mp.Aligner(index)  
inputFile = open(inputFile, 'r')

## Create output bed file with header
outPath = outDir + '/transducedDb.bed'
outFile = open(outPath, 'w') 

row = 'chrom' + "\t" + 'tdBeg' + "\t" + 'tdEnd' + "\t" + 'cytobandId' + '\n'
outFile.write(row)

## For each source element
for line in inputFile:
	line = line.rstrip('\r\n')
    
	## Discard header
	if not line.startswith("#"):        
		fieldsList = line.split("\t")
		chrom, beg, end, name, cytobandId, score, strand = fieldsList
		
		## a) Element in plus
		# ---------------> end ........transduced........ end + 15000
		if (strand == '+'):

			# retrieve a subsequence from the index
			tdBeg = int(end)
			tdEnd = int(end) + 15000

			sequence = index.seq(chrom, tdBeg, tdEnd)     

		## b) Element in minus
		# beg - 15000 ........transduced........ beg <--------------- end  
		else:

			# retrieve a subsequence from the index
			tdBeg = int(beg) - 15000
			tdEnd = int(beg)

			sequence = index.seq(chrom, tdBeg, tdEnd) 

		coord = chrom + ':' + str(tdBeg) + '-' + str(tdEnd)

		## Add source element transduced region to the fasta
		header = 'transduced|L1|' + cytobandId + '|' + coord
		FASTA.seqDict[header] = sequence

		## Write transduced region coordinates into bed file
		row = chrom + "\t" + str(tdBeg) + "\t" + str(tdEnd) + "\t" + cytobandId + '\n'
		outFile.write(row)
        
## 2. Write fasta file
########################
filePath = outDir + '/transducedDb.fa'
FASTA.write(filePath)

## 3. Mask interspersed repeats in the transduced regions sequence
####################################################################

## Create logs directory
logDir = outDir + '/Logs'
unix.mkdir(logDir)

## Do masking
err = open(logDir + '/RepeatMasker.err', 'w') 
out = open(logDir + '/RepeatMasker.out', 'w') 
command = 'RepeatMasker -nolow -norna -no_is -div 18 -species human ' + filePath 
status = subprocess.call(command, stderr=err, stdout=out, shell=True)


if status != 0:
	step = 'REPEATMASKER'
	msg = 'Repeat masking failed' 
	log.step(step, msg)

print('***** Finished! *****')
print()