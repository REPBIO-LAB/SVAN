'''
Module 'assembly' - Contains funtions to create contigs from a set of input sequences
'''

## DEPENDENCIES ##
# External
import subprocess
import os

# Internal
from GAPI import log
from GAPI import unix
from GAPI import alignment
from GAPI import formats


## FUNCTIONS ##
    
def assemble_overlap(fastaA, fastaB, technology, outDir):
    '''
    Assemble two sets sequences via the identification of reciprocal overlap at the sequence ends. 
    (custom algorithm based on minimap2 alignments)

     -: aligned; _: clipped      sequence_bkp
    sequence A:     -------------------*_____________________
    sequence B:                                    |________|______________*-------------------
    assembled:      ---------------------------------------------------------------------------

    Input:
        1. fastaA: FASTA file containing set of sequences A
        2. fastaB: FASTA file containing set of sequences B
        3. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        4. outDir: output directory

    Output:
        1. contigFile: FASTA file containing the generated contig or 'None' if no reciprocal overlap is found 
    '''    
    ## 0. Create directories ##
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Align sequences A vs B ##
    ## Set preset according to the technology
    preset = alignment.minimap2_presets(technology)
    
    ## Do alignment 
    PAF_file = outDir + '/overlaps.paf'
    err = open(logDir + '/minimap2.err', 'w') 
    command = 'minimap2 -x ' + preset + ' ' + fastaB + ' ' + fastaA + ' > ' + PAF_file

    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ASSEMBLE-OVERLAP'
        msg = 'Sequences A vs B alignment failed' 
        log.step(step, msg)

    ## 2. Read PAF alignments ##
    PAF = formats.PAF()
    PAF.read(PAF_file)

    # Exit function if no hit found
    if not PAF.alignments:
        return None

    ## 3. Filter overlaps 
    # Filtering criteria:
    # Pick alignments on the + strand 
    # Pick alignments starting within 250 from sequences ends. 
    filtered = []

    for overlap in PAF.alignments:

        # Compute the distance between the corresponding sequence end and the overlap begin or end position
        distA = overlap.qLen - overlap.qEnd 
        distB = overlap.tBeg
    
        # Apply filters
        if (overlap.strand == '+') and (distA <= 250) and (distB <= 250):
            filtered.append(overlap)

    # Replace raw by filtered overlaps
    PAF.alignments = filtered

    # Exit function if all the overlapping hits have been filtered
    if not PAF.alignments:
        return None

    ## 4. Pick longest overlap passing the filters
    # Sort PAF in decreasing overlap lengths
    PAF.alignments = PAF.sortByLen()

    # Pick longest overlap
    overlapLongest = PAF.alignments[0]

    ## 5. Concatenate overlapping sequences to generate a contig 
    ## Read input fasta files
    sequencesA = formats.FASTA()
    sequencesA.read(fastaA)
    
    sequencesB = formats.FASTA()
    sequencesB.read(fastaB)
    
    ## Contig generation  
    #  -: aligned; _: clipped                     sequence_bkp    qBeg       qEnd
    # sequence A (query):                ---------------*___________|_________|___      sequence_bkp
    # sequence B (template):                                    ____|_________|______________*-------------------
    #                                                             tBeg       tEnd
    # contig:                            -------------------------------------|----------------------------------
    #                                                  sequence_1        :qEnd tEnd:       sequence_2
    seqA = sequencesA.seqDict[overlapLongest.qName][:overlapLongest.qEnd]
    seqB = sequencesB.seqDict[overlapLongest.tName][overlapLongest.tEnd:]
    contigSeq = seqA + seqB 

    ## Write contig into fasta
    contig = formats.FASTA()
    contig.seqDict['CONTIG'] = contigSeq

    contigFile = outDir + '/contig.fa'
    contig.write(contigFile)

    return contigFile

def polish_racon(templates, sequences, technology, nbRounds, outDir):
    '''
    Use a collection of sequences to polish a set of target sequences (i.e. assembled contig) 
    
    Input:
        1. templates: FASTA file containing sequences to be polished
        2. sequences: FASTA file containing set of sequences used to polish the templates
        3. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        4. nbRounds: number of polishing rounds to be performed 
        5. outDir: output directory
        
    Output:
        1. polished: FASTA file containing polished sequences or 'None' if pipeline fails at any step
    '''
    ## Create logs directory:
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## Set preset according to the technology
    preset = alignment.minimap2_presets(technology)

    ## Start with polished set as None
    polished = None            

    ## For each polishing round
    for roundId in range(1, nbRounds + 1):
        
        ## 1. Align reads against the template ##
        PAF = outDir + '/alignments_' + str(roundId) + '.paf'
        err = open(logDir + '/minimap2_' + str(roundId) + '.err', 'w') 
        command = 'minimap2 -x ' + preset + ' ' + templates + ' ' + sequences + ' > ' + PAF
        status = subprocess.call(command, stderr=err, shell=True)

        if status != 0:
            step = 'POLISH-RACON'
            msg = 'Alignment of sequences against template failed' 
            log.step(step, msg)

            polished = None            
            break

        ## 2. Template polishing with racon ##
        polished = outDir + '/polished_' + str(roundId) + '.fa'
        err = open(logDir + '/racon_' + str(roundId) + '.err', 'w') 
        command = 'racon --include-unpolished ' + sequences + ' ' + PAF + ' ' + templates + ' > ' + polished
        status = subprocess.call(command, stderr=err, shell=True)

        if status != 0:
            step = 'POLISH-RACON'
            msg = 'Template polishing failed' 
            log.step(step, msg)

            polished = None            
            break

        ## 3. Set polished as templates prior attempting a new polishing round
        templates = polished

    return polished


def getConsensusSeq(FASTA_file, outDir):
    '''
    Build consensus seq from fasta file

    Input:
        1. FASTA_file: Fasta file
        2. outDir
    Output:
        1. consensusPath: consensus file
        2. consensusSeq: consensus sequence
    '''

    # 1. Check that the fasta file is not empty:
    if not os.stat(FASTA_file).st_size == 0:
        
        # 2. Make multiple sequence alignment
        msfPath = FASTA_file.replace("fa", "msf")
        # TODO: Add muscle to environment
        err = open(outDir + '/muscle.err', 'w') 
        command = 'muscle -in ' + FASTA_file + ' -out ' + msfPath + ' -msf'

        status = subprocess.call(command, stderr=err, shell=True)

        if status != 0:
            step = 'MUSCLE'
            msg = 'Muscle alignment failed' 
            log.step(step, msg)

        # 3. Generate consensus sequence (cons tool from EMBOSS package)
        consensusPath = FASTA_file.replace("_supportingReads", "_consensus")
        # TODO: Add emboss to environment
        err = open(outDir + '/cons.err', 'w') 
        command = 'cons -sequence ' + msfPath + ' -outseq ' + consensusPath + ' -identity 0 -plurality 0'
        status = subprocess.call(command, stderr=err, shell=True)

        if status != 0:
            step = 'CONS'
            msg = 'Cons failed' 
            log.step(step, msg)

        # 4. Check that consensus file exists and is not empty
        if os.path.exists(consensusPath) and not os.stat(consensusPath).st_size == 0:
            
            # 5. Read consensus sequence 
            consensusFastaObj = formats.FASTA()
            consensusFastaObj.read(consensusPath)
            consensusSeq = consensusFastaObj.seqDict["EMBOSS_001"].upper()

            # 6. Replace '-' by 'N' for ambiguous bases:
            consensusSeq = consensusSeq.replace('-', 'N')

            # 7. Convert consensus sequence into upper case:
            consensusSeq = consensusSeq.upper()
            
        else:
            consensusSeq = None

    return consensusPath, consensusSeq

