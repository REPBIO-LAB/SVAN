#!/usr/bin/env python
#coding: utf-8

## DEPENDENCIES ##
# External
import argparse
import os
import sys

# Internal
from GAPI import formats


######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('repeatMaskerDb', help='')
parser.add_argument('subfamilies', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='Output directory. Default: current working directory')

## 2. Parse user´s input and initialize variables ##
args = parser.parse_args()
repeatMaskerDb = args.repeatMaskerDb
subfamilies = args.subfamilies
outDir = args.outDir


##############################################
## Display configuration to standard output ##
##############################################
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('repeatMaskerDb: ', repeatMaskerDb)
print('subfamilies: ', subfamilies)
print('outDir: ', outDir, "\n")

print('***** Executing ', scriptName, '.... *****', "\n")


##########
## CORE ## 
##########

## 1. Read repeat masker database
##################################

## Initialize FASTA object
FASTA = formats.FASTA()

## Read fasta
FASTA.read(repeatMaskerDb)

## 2. Create list of target subfamilies
########################################
subfamilies = open(subfamilies, 'r')

targets = {}

## For each target subfamily
for line in subfamilies:
    line = line.rstrip('\r\n')
    
    ## Discard header
    if not line.startswith("#"):        
        fieldsList = line.split("\t")
        subfamily, family = fieldsList
        targets[subfamily] = family

## 3. Select target subfamilies consensus sequences
####################################################
FASTA_out = formats.FASTA()

# For each consensus sequence in the repeatmasker database
for sequenceId in FASTA.seqDict:
    subfamily = sequenceId.split('#')[0]

    # Consensus sequence is a target
    if subfamily in targets:

        family = targets[subfamily]
        newId = 'consensus|'+ family + '|' + subfamily  
        sequence = FASTA.seqDict[sequenceId]

        ## Remove trailing polyA from consensus sequence (Include polyA trimming as option)
        nbTrimmedBp = 0

        for nucleotide in reversed(sequence):
            # Trailing A
            if (nucleotide == 'A') or (nucleotide == 'a'):
            	nbTrimmedBp += 1
            # Not A, stop trimming
            else:
            	break
        
        index = len(sequence) - nbTrimmedBp
        trimmed = sequence[:index]
        FASTA_out.seqDict[newId] = trimmed

        #FASTA_out.seqDict[newId] = sequence



## 4. Write fasta file
########################
filePath = outDir + '/consensusDb.fa'
FASTA_out.write(filePath)

print('***** Finished! *****')
print()
