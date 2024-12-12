'''
Module 'alignment' - Contains funtions align sequences
'''

## DEPENDENCIES ##
# External
import subprocess
import os
import pysam

# Internal
from GAPI import log
from GAPI import unix
from GAPI import bamtools 
from GAPI import formats


## FUNCTIONS ##
def sam2paf(SAM, fileName, outDir):
    '''
    Wrapper to convert a SAM into a PAF file
    '''
    PAF = outDir + '/' + fileName + '.paf'
    err = open(outDir + '/sam2paf.err', 'w') 
    command = 'paftools.js sam2paf -p ' + SAM + ' > ' + PAF
    subprocess.call(command, stderr=err, shell=True)

    return PAF    


def index_minimap2(fastaPath, fileName, outDir):
    '''
    Wrapper to generate minimap2 index for fasta file
    '''
    indexPath = outDir + '/' + fileName + '.mmi'
    err = open(outDir + '/index.err', 'w') 
    command = 'minimap2 -k 11 -w 1 -d ' + indexPath + ' ' + fastaPath 
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'INDEX'
        msg = 'minimap2 indexing failed' 
        log.step(step, msg)
    
    return indexPath

def index_bwa(fastaPath, outDir):
    '''
    Wrapper to generate BWA index for fasta file
    '''
    err = open(outDir + '/index.err', 'w') 
    command = 'bwa index ' + fastaPath 
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'INDEX'
        msg = 'BWA indexing failed' 
        log.step(step, msg)
    
def minimap2_presets(technology):
    '''
    Set minimap2 preset according to the sequencing technology

    Input:
        1. technology: Sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        
    Output:
        1. preset: proper preset  
    '''

    if technology == 'PACBIO':
        preset = 'map-pb'
        
    elif technology == 'NANOPORE':
        preset = 'map-ont'

    elif technology == 'ILLUMINA':
        preset = 'sr'

    else:
        log.info('ERROR: ' + technology + ' technology not supported')
        sys.exit(1)

    return preset

def alignment_minimap2(FASTA, index, fileName, processes, outDir):
    '''
    Align a set of sequence into a reference with minimap2

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. index: Path to the the index of the reference in .mmi format (generated with minimap2)
        3. fileName: output file will be named accordingly
        4. processes: Number of processes used by minimap2
        5. outDir: Output directory

    Output:
        1. PAF: Path to PAF file containing input sequences alignments or 'None' if alignment failed 
    '''
    ## Align the sequences into the reference
    # Note, condider to use -Y to get soft clippings for supplementary alignments
    PAF = outDir + '/' + fileName + '.paf'
    err = open(outDir + '/align.err', 'w') 
    command = 'minimap2 -t ' + str(processes) + ' ' + index + ' ' + FASTA + ' > ' + PAF
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Local alignment failed' 
        log.step(step, msg)

    return PAF

def alignment_minimap2_hifi(FASTA, ref, fileName, processes, outDir):
    '''
    Align a set of sequence into a reference with minimap2

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. ref: Path to the fasta of the reference sequence 
        3. fileName: output file will be named accordingly
        4. processes: Number of processes used by minimap2
        5. outDir: Output directory

    Output:
        1. PAF: Path to PAF file containing input sequences alignments or 'None' if alignment failed 
    '''
    ## Align the sequences into the reference
    # Note, condider to use -Y to get soft clippings for supplementary alignments
    PAF = outDir + '/' + fileName + '.paf'
    err = open(outDir + '/align.err', 'w') 
    command = 'minimap2 -x map-hifi -t ' + str(processes) + ' ' + ref + ' ' + FASTA + ' > ' + PAF
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Local alignment failed' 
        log.step(step, msg)

    return PAF


def alignment_minimap2_bam(FASTA, index, fileName, processes, outDir):
    '''
    Align a set of sequence into a reference with minimap2

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. index: Path to the the index of the reference in .mmi format (generated with minimap2)
        3. fileName: output file will be named accordingly
        4. processes: Number of processes used by minimap2
        5. outDir: Output directory

    Output:
        1. BA;: Path to BAM file containing input sequences alignments or 'None' if alignment failed 
    '''
    ## 1. Align the sequences into the reference
    # Note, condider to use -Y to get soft clippings for supplementary alignments
    SAM = outDir + '/' + fileName + '.sam'
    err = open(outDir + '/align.err', 'w') 
    command = 'minimap2 -a -t ' + str(processes) + ' ' + index + ' ' + FASTA + ' > ' + SAM
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Local alignment failed' 
        log.step(step, msg)

    ## 2. SAM to BAM conversion
    BAM = bamtools.SAM2BAM(SAM, outDir)

    return BAM


def alignment_minimap2_spliced(FASTA, index, fileName, processes, outDir):
    '''
    Align a set of sequence into a reference with minimap2

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. index: Path to the the index of the reference in .mmi format (generated with minimap2)
        3. fileName: output file will be named accordingly
        4. processes: Number of processes used by minimap2
        5. outDir: Output directory

    Output:
        1. SAM: Path to SAM file containing input sequences alignments or 'None' if alignment failed 
    '''
    ## Align the sequences into the reference
    # Note, condider to use -Y to get soft clippings for supplementary alignments
    SAM = outDir + '/' + fileName + '.sam'
    err = open(outDir + '/align.err', 'w') 
    command = 'minimap2 -a -x splice -t ' + str(processes) + ' ' + index + ' ' + FASTA + ' > ' + SAM
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Local alignment failed' 
        log.step(step, msg)

    return SAM

def alignment_bwa(FASTA, reference, fileName, processes, outDir):
    '''
    Align a set of sequence into a reference with bwa mem

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. reference: Path to the reference genome in fasta format (bwa mem index must be located in the same folder)
        3. fileName: output file will be named accordingly
        4. processes: Number of processes used by bwa mem(
        5. outDir: Output directory

    Output:
        1. SAM: Path to SAM file containing input sequences alignments or 'None' if alignment failed 
    '''
    ## Align the sequences into the reference
    SAM = outDir + '/' + fileName + '.sam'
    err = open(outDir + '/align.err', 'w') 
    command = 'bwa mem -k 8 -T 0 -Y -t ' + str(processes) + ' ' + reference + ' ' + FASTA + ' > ' + SAM
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Alignment failed' 
        log.step(step, msg)

    return SAM

def alignment_bwa_2(FASTA, reference, fileName, processes, outDir):
    '''
    Align a set of sequence into a reference with bwa mem

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. reference: Path to the reference genome in fasta format (bwa mem index must be located in the same folder)
        3. fileName: output file will be named accordingly
        4. processes: Number of processes used by bwa mem
        5. outDir: Output directory

    Output:
        1. SAM: Path to SAM file containing input sequences alignments or 'None' if alignment failed 
    '''
    ## Align the sequences into the reference
    SAM = outDir + '/' + fileName + '.sam'
    err = open(outDir + '/align.err', 'w') 
    command = 'bwa mem -Y -t ' + str(processes) + ' ' + reference + ' ' + FASTA + ' > ' + SAM
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Alignment failed' 
        log.step(step, msg)

    return SAM

# NOTE MERGE SR2020: OLD IN MASTER
def alignment_blat_oldMaster(FASTA, reference, fileName, outDir):
    '''
    Align a set of sequence into a reference with blat

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. reference: Path to the reference genome in fasta format 
        3. fileName: output file will be named accordingly
        4. outDir: Output directory

    Output:
        1. SAM: Path to SAM file containing input sequences alignments or 'None' if alignment failed 
    '''
    ## Align the sequences into the reference
    PSL = outDir + '/' + fileName + '.psl'
    err = open(outDir + '/align.err', 'w') 
    command = 'blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -noHead -out=psl ' + reference + ' ' + FASTA + ' ' + PSL
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Alignment failed' 
        log.step(step, msg)

    return PSL

def alignment_blat(FASTA, reference, args, fileName, outDir):
    '''    
    Align a set of sequence into a reference with blat

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. reference: Path to the reference genome in fasta format (bwa mem index must be located in the same folder)
        3. args: dictionary containing blat arguments
        4. fileName: output file will be named accordingly
        5. outDir: Output directory

    Output:
        1. SAM: Path to SAM file containing input sequences alignments or 'None' if alignment failed 
    '''
    ## Align the sequences into the reference
    PSL = outDir + '/' + fileName + '.psl'
    err = open(outDir + '/align.err', 'w') 

    # Set blat arguments
    blatArgs = []

    if 'stepSize' in args.keys():
        blatArgs.append('-stepSize='+str(args['stepSize']))

    if 'tileSize' in args.keys():
        blatArgs.append('-tileSize='+str(args['tileSize']))

    command = 'blat ' + ' '.join(blatArgs) + ' -repMatch=2253 -minScore=20 -minIdentity=0 -noHead -out=psl ' + reference + ' ' + FASTA + ' ' + PSL
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Alignment failed' 
        log.step(step, msg)

    return PSL

def targeted_alignment_minimap2(FASTA, targetInterval, reference, outDir, outFormat):
# NOTE 2020: In 2020:
# def targeted_alignment_minimap2(FASTA, targetInterval, reference, outDir):
    '''
    Align a set of sequences into a reference target region. 
    
    Useful for doing local realignment of reads around SV breakpoints. Much faster than whole genome realignment

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. targetInterval: Reference genome interval where sequences will be aligned. The interval must be provided as chr:beg-end.
        3. reference: Path to the reference sequences in fasta format. An index of the reference generated with samtools faidx must be located in the same directory
        4. outDir: Output directory
        5. outFormat: BAM or SAM

    Output:
        1. BAM: Path to sorted BAM file containing input sequences alignments or 'None' if realignment failed 
    '''
    ## 0. Create logs directory
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Extract the reference target region prior alignment 
    target = outDir + '/target.fa'
    err = open(logDir + '/target.err', 'w') 
    command = 'samtools faidx ' + reference + ' ' + targetInterval + ' > ' + target
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'TARGET'
        msg = 'Extraction of reference target region failed' 
        log.step(step, msg)
        return None

    ## 2. Align the sequences into the target region 
    # Use -Y to get soft clippings for supplementary alignments
    SAM = outDir + '/alignments.sam'
    err = open(logDir + '/align.err', 'w') 
    command = 'minimap2 -Y -a ' + target + ' ' + FASTA + ' > ' + SAM
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Local alignment failed' 
        log.step(step, msg)
        return None

    # NOTE 2020: In 2020 these 2 lines were deleted:
    if outFormat == "SAM":
        return SAM

    ## 3. Convert SAM to sorted BAM
    BAM = bamtools.SAM2BAM(SAM, outDir)

    ## 4. Do cleanup
    unix.rm([target, SAM])

    return BAM    

def targetered2genomic_coord(event, ref, offset):
    '''
    Convert event coordinates resulting from the realignment of a sequence into a target region into genomic coordinates

    Input:
        1. event: INS, DEL or CLIPPING event object
        2. ref: reference corresponding to the targetered seq
        3. offset: offset to be added to event coordinates

    Output:
        1. event: modified event object
    '''
    ## 1. Convert reference                
    event.ref = ref

    ## 2. Add offset to event begin and end coordinates
    event.beg = event.beg + offset
    event.end = event.end + offset 

    ## 3. Add offset to supplementary alignments 
    if event.supplAlignment != None:

        supplAlignments = ''

        # For each supplementary alignment
        for supplAlignment in event.supplAlignment.split(';'):

            # Stop after the last supplementary alignment
            if supplAlignment == '':
                break

            supplRef, supplBeg, strand, cigar, MAPQ, NM = supplAlignment.split(',')

            supplRef = ref
            beg = int(supplBeg) + offset

            supplAlignment = supplRef + ',' + str(beg) + ',' + strand + ',' + cigar + ',' + MAPQ + ',' + NM
            supplAlignments = supplAlignments + supplAlignment + ';' if supplAlignments != '' else supplAlignment + ';'

            event.supplAlignment = supplAlignments

    return event

def organize_hits_paf(PAF_path):
    '''
    Group hits by query name into a dictionary

    Input:
        1. PAF_path: Path to paf file containing alignments

    Output:
        1. hits: dictionary containing query names as keys and the list of alignments for each query as values
    '''
    ## 1. Read PAF 
    PAF = formats.PAF()
    PAF.read(PAF_path)

    ## 2. Organize hits by query name into the dictionary
    hits = {}

    # For each read alignment 
    for alignment in PAF.alignments:

        # First hit for this query, initialize PAF
        if alignment.qName not in hits:
            PAF = formats.PAF()
            hits[alignment.qName] = PAF

        # Add hit to PAF
        hits[alignment.qName].alignments.append(alignment)

    return hits     

def getPAFAlign(FASTA_file, indexDb, outDir):
    '''
    Align FASTA with minimap2 and retrieve PAF file
    Input:
        1. FASTA_file
        2. indexDb
        3. outDir
    Output:
        1. PAF file of alignment
    '''
    # Name of output PAF file
    PAF_file = FASTA_file.replace(".fa", "_alignments.paf")

    # Align FASTA against given db 
    command = 'minimap2 ' + indexDb + ' ' + FASTA_file + ' > ' + PAF_file
    err = open(outDir + '/minimap2.err', 'w') 
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN-INSERT'
        msg = 'minimap2 alignment failed' 
        log.step(step, msg)
    
    return PAF_file