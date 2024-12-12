'''
Module 'repeats' - Contains functions for the identification and characterization of simple repeats from sequences
'''

## DEPENDENCIES ##
# External
import subprocess

# Internal
from GAPI import log
from GAPI import unix
from GAPI import formats
from GAPI import alignment
from GAPI import sequences

## FUNCTIONS ##
def is_simple_repeat(FASTA_file, minPercSimple, outDir):
    '''
    Determine if an input sequence is a simple repeat
    
    Input:
        1. FASTA_file: Path to FASTA file containing the sequence
        2. minPercSimple: Minimum percentage of bases corresponding to simple repeats for classifying the sequence as simple repeat 
        3. outDir: Output directory
        
    Output:
    ''' 

    ## 0. Create logs directory 
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Run dustmasker on the input sequence 
    output = outDir + '/dustmasker.out'
    command = 'dustmasker -outfmt acclist -in ' + FASTA_file + ' -out ' + output
    status = subprocess.call(command, shell=True)

    if status != 0:
        step = 'DUSTMASKER'
        msg = 'Dustmasker failed' 
        log.step(step, msg)

    ## 2. Compute the total number of bp corresponding to simple repeats
    output = open(output)
    simpleLen = 0  

    # For line in the file
    for line in output:
            
        # Select lines starting with '>' 
        if line.startswith('>'):
            line = line.rstrip()

            ## Compute the length of the simple repeat hit
            beg, end = line.split('\t')[1:]
            hitLen = int(end) - int(beg)

            ## Add the hit length to the total simple repeat length
            simpleLen += hitLen
        
    ## 3. Compute the percentage of bp corresponding to simple repeat DNA
    ## First compute input sequence length
    FASTA = formats.FASTA()
    FASTA.read(FASTA_file)
    sequence = list(FASTA.seqDict.values())[0]
    seqLen = len(sequence)

    ## Then the percentage
    percSimple = float(simpleLen) / seqLen * 100
    
    ## 4. Classify sequence as simple repeat if at least X% of its sequence correspond to a simple repeat 
    if (percSimple >= minPercSimple):
        status = 'resolved'
        insType = 'simple_repeat'

    else:
        status = 'unresolved'        
        insType = None

    return insType, status, percSimple
