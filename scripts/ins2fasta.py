## DEPENDENCIES ##
# External
import os
import sys
import argparse

# Internal
import formats
import sequences
import retrotransposons

######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('vcf', help='')
parser.add_argument('outDir', help='Output directory')

## 2. Parse user input ##
args = parser.parse_args()
vcf = args.vcf
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('vcf: ', vcf)
print('outDir: ', outDir, "\n")


##########
## CORE ##
##########

## 1. Load vcf file
####################
VCF = formats.VCF()
VCF.read(vcf)

## 2. Create fasta object
#########################

## Initialize fasta
fasta = formats.FASTA()

## Add sequences for full L1s to fasta
for variant in VCF.variants:

    insId = variant.ID + '_' + variant.chrom + ':' + str(variant.pos) 

    fasta.seqDict[insId] = variant.alt

## 3. Write fasta file as output
################################
outFile = outDir + '/insertions_seq.fa'
fasta.write(outFile)
