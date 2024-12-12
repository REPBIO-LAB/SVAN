'''
Module 'bamtools' - Contains functions for extracting data from bam files
'''

## DEPENDENCIES ##
# External
import pysam
import subprocess
import sys
from cigar import Cigar
import numpy as np
import statistics

# Internal
from GAPI import log
from GAPI import unix
from GAPI import events
from GAPI import gRanges
from GAPI import formats
from GAPI import sequences

## FUNCTIONS ##
def get_refs(bam):
    '''
    Get all references present in the bam file.

	Input:
		1. bam: indexed BAM file
	
	Output:
		1. refs: String containing all references from the bam file, separated by commas.
    '''
    bamFile = pysam.AlignmentFile(bam, 'rb')
    refs  = ','.join(bamFile.references)
    return refs


def get_ref_lengths(bam):
    '''
    Make dictionary containing the length for each reference

	Input:
		1. bam: indexed BAM file
	
	Output:
		1. lengths: Dictionary containing reference ids as keys and as values the length for each reference
    '''
    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, 'rb')
    
    ## Make dictionary with the length for each reference (move this code into a function)
    refLengths = dict(list(zip(bamFile.references, bamFile.lengths)))
    
    ## Close bam file
    bamFile.close()
    
    return refLengths

def infer_readSize(bam):
    '''
    Infer read size from bam file

    Note: we should improve the script by taking 1000 random reads from the bam file
    '''
        
    # take the first 500 reads of the bam file
    command = 'samtools view ' + bam + '| awk \'{print length($10)}\' | head -500 | tr \'\n\' \' \''
    result = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        
    # if command fails, exit
    if result.returncode != 0:
        step = 'infer_readSize'
        msg = 'readSize inference failed' 
        log.step(step, msg)
        sys.exit(1)
        
    # save the result in a list of integers
    readSizes_str = result.stdout.decode('utf-8').split(" ")
    readSizes_str.remove("")
    readSizes_int = [int(i) for i in readSizes_str]
        
    # calculate the mode
    readSize = statistics.mode(readSizes_int)
        
    return(readSize)

def alignment_length_cigar(CIGAR):
    '''
    Compute alignment on the reference length from CIGAR string

    Input:
        1. CIGAR: CIGAR string

    Output:
        1. alignmentLen: alignment on the reference length
    '''
    ## 1. Read CIGAR string using proper module
    cigarTuples = Cigar(CIGAR)

    ## 2. Iterate over the operations and compute the alignment length
    alignmentLen = 0

    for cigarTuple in list(cigarTuples.items()):

        length = int(cigarTuple[0])
        operation = cigarTuple[1]

        ### Update reference alignment length
        ## a) Operations consuming query and reference
        # - Op M, tag 0, alignment match (can be a sequence match or mismatch)
        # - Op =, tag 7, sequence match
        # - Op X, tag 8, sequence mismatch
        if (operation == 'M') or (operation == '=') or (operation == 'X'):
            alignmentLen += length

        ## b) Operations only consuming reference
        # - Op D, tag 2, deletion from the reference
        # - Op N, tag 3, skipped region from the reference
        elif (operation == 'D') or (operation == 'N'):
            alignmentLen += length
            
    return alignmentLen


def alignment_interval_query(CIGAR, orientation):
    '''
    Compute alignment on the reference length from CIGAR string

    Input:
        1. CIGAR: CIGAR string
        2. orientation: alignment orientation (+ or -) 

    Output:
        1. beg: begin position in query
        2. end: end position in query
    '''
    ## 1. Read CIGAR string using proper module
    cigar = Cigar(CIGAR)

    ## 2. Iterate over the operations and compute query alignment length and start position in query
    alignmentLen = 0
    counter = 0 # Count operations

    for cigarTuple in list(cigar.items()):

        length = int(cigarTuple[0])
        operation = cigarTuple[1]

        ## Set start position in query based on first operation 
        if counter == 0:

            # a) Soft or Hard clipping
            if (operation == 'S') or (operation == 'H'):
                startPos = length

            # b) No clipping
            else:
                startPos = 0
            
        #### Update query alignment length
        # - Op M, alignment match (can be a sequence match or mismatch)
        # - Op =, sequence match
        # - Op X, sequence mismatch
        # - Op I, insertion to the reference
        if (operation == 'M') or (operation == '=') or (operation == 'X') or (operation == 'I'):
            alignmentLen += length

        ## Update operations counter
        counter += 1

    ## 3. Compute alignment interval in raw query
    ## Compute read length
    readLen = len(cigar)

    # a) Query aligned in +
    if orientation == '+':
        beg = startPos
        end = startPos + alignmentLen

    # b) Query aligned in - (reversed complemented to align)
    else:
        beg = readLen - startPos - alignmentLen
        end = readLen - startPos
        
    return beg, end


def SAM2BAM(SAM, outDir):
    '''
    Convert SAM file into sorted BAM and make BAM index

	Input:
		1. SAM: File containing alignments in SAM format

	Output:
		1. BAM_sorted: Sorted and indexed BAM file. BAM index located in the same directory with the extension '.bai'
    '''
    ## 0. Create logs directory 
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Convert SAM into BAM 
    BAM = outDir + '/alignments.bam'
    err = open(logDir + '/SAM2BAM.err', 'w') 
    command = 'samtools view -Sb ' + SAM + ' > ' + BAM
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'SAM2BAM'
        msg = 'SAM to BAM conversion failed' 
        log.step(step, msg)

    ## 2. Sort bam 
    BAM_sorted = outDir + '/alignments.sorted.bam'
    err = open(logDir + '/sort.err', 'w') 
    command = 'samtools sort ' + BAM + ' > ' + BAM_sorted
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'SORT'
        msg = 'BAM sorting failed' 
        log.step(step, msg)    

    ## 3. Index bam 
    BAM_index = outDir + '/alignments.sorted.bam.bai'
    err = open(logDir + '/index.err', 'w') 
    command = 'samtools index ' + BAM_sorted + ' > ' + BAM_index
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'INDEX'
        msg = 'BAM indexing failed' 
        log.step(step, msg)

    return BAM_sorted

def BAM2BED(BAM, outDir):
    '''
    Convert BAM file into BED using bedtools

	Input:
		1. BAM: Path to BAM file 
        2. outDir: Output directory

	Output:
		1. BED: Path to BED file
    '''
    ## 0. Create logs directory 
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Convert BAM into BED
    BED_path = outDir + '/alignments.bed'
    err = open(logDir + '/BAM2BED.err', 'w') 
    command = 'bedtools bamtobed -split -i ' + BAM + ' > ' + BED_path
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'BAM2BED'
        msg = 'BAM to BED conversion failed' 
        log.step(step, msg)
    
    ## 2. Add header to BED file
    header = "#ref \t beg \t end \t name \t score \t strand \n"
    with open(BED_path, 'r') as original: data = original.read()
    with open(BED_path, 'w') as modified: modified.write(header + data)

    return BED_path

def phred2ASCII(phred):
    '''
    Convert Phred quality scores into ASCII (Sanger format used in FASTQ)

	Input:
		1. phred: List containing per base standard Phred quality scores (from 0 to 93) as provided by pysam.AlignmentFile.query_qualities attribute.
	
	Output:
		1. ASCII: List containing per base phred quality scores encoded by ASCII characters 33 to 126 
    '''
    ASCII = [chr(x + 33) for x in phred]

    return ASCII


def BAM2FASTQ_entry(alignmentObj):
    '''
    Transform a BAM alignment into a FASTQ_entry object. 

	Input:
		1. alignmentObj: pysam.AlignedSegment object.
	
	Output:
		1. FASTQ_entry: formats.FASTQ_entry object
    '''
    ## 1. Convert Phred quality scores to ASCII 
    if not alignmentObj.query_qualities == None:
        ASCII = phred2ASCII(alignmentObj.query_qualities)
        ASCII = "".join(ASCII)
    else:
        qual = None

    ## 2. Obtain raw read and quality strings (Prior alignment)
    # a) Read mapped in reverse -> Make complementary reverse of the sequence and the reverse of the quality 
    if alignmentObj.is_reverse:
        seq = sequences.rev_complement(alignmentObj.query_sequence)

        if not alignmentObj.query_qualities == None:
            qual = ASCII[::-1]
                
    # b) Read mapped in forward
    else:
        seq = alignmentObj.query_sequence

        if not alignmentObj.query_qualities == None:
            qual = ASCII

    ## 3. Create FASTQ_entry object
    FASTQ_entry = formats.FASTQ_entry(alignmentObj.query_name, seq, '', qual)
    return FASTQ_entry

def alignments2PAF(alignments):
    '''
    Convert as set of pysam aligned segments into a PAF object

	Input:
		1. alignments: list of aligned segments

	Output:
		1. PAF: PAF object containing alignments
    '''
    
    ## 1. Initialize PAF object
    PAF = formats.PAF()

    ## 2. Convert each aligned segment into a PAF_alignment object and add to PAF
    for alignment in alignments:

        # Discard unmapped sequences
        if not alignment.is_unmapped:
            
            strand = '-' if alignment.is_reverse else '+'
            fields = [alignment.query_name, alignment.infer_read_length(), alignment.query_alignment_start, alignment.query_alignment_end, strand, alignment.reference_name, alignment.reference_length, alignment.reference_start, alignment.reference_end, 0, 0, alignment.mapping_quality]
            line = formats.PAF_alignment(fields)
            PAF.alignments.append(line)

    return PAF


def binning(targetBins, bam, binSize, targetRefs):
    '''
    Split the genome into a set of genomic bins. Two possible binning approaches:
    1) Use predefined bins if bed file provided (targetBins)
    2) Non overlapping bins of a given size (binSize) for a set of target references (targetRefs). Reference length extracted from the input 'bam' file

    Input:
        1. targetBins: Bed file containing predefined bins OR None (in this case bins will be created de novo)
        2. bam: BAM file used to know the length of the target references 
        3. binSize: Binning size
        3. targetRefs: Comma separated list of target references

    Output:
        1. bins: List of bins. Each list item corresponds to a list [ref, beg, end]
    '''

    # A) Create bins de novo
    if targetBins == None:

        ## Split the reference genome into a set of genomic bins
        bins = makeGenomicBins(bam, binSize, targetRefs)

    # B) Read bins from bed file
    else:
        BED = formats.BED()
        BED.read(targetBins, 'List', None)        
        bins = [ [line.ref, line.beg, line.end] for line in BED.lines]
    
    return bins


def makeGenomicBins(bam, binSize, targetRefs):
    '''
    Split the genome into a set of non overlapping bins of 'binSize' bp.

    Input:
        1. bam: BAM file
        2. binSize: size of the bins
        3. targetRefs: list of target references. None if all the references are considered

    Output:
        1. bins: List of non overlapping bins. Each list item corresponds to a tuple (ref, beg, end)
    '''
    ## Obtain the length for each reference
    refLengths = get_ref_lengths(bam)

    ## Select target references
    if targetRefs != None:
        targetRefs = [str(i) for i in targetRefs]
        refLengths  = {ref: refLengths [ref] for ref in targetRefs}

    ## Split each reference into evenly sized bins
    bins = []

    # For each reference
    for ref, length in refLengths .items():

        ## Define bins boundaries
        boundaries = [boundary for boundary in range(0, length, binSize)]
        boundaries = boundaries + [length]

        ## Make bins
        for idx, beg in enumerate(boundaries):

            ## Skip last element from the list
            if beg < boundaries[-1]:
                end = boundaries[idx + 1]
                window = (ref, beg, end)
                bins.append(window)

    return bins


def collectSV_paired(ref, binBeg, binEnd, tumourBam, normalBam, confDict):
    '''
    Collect structural variant (SV) events in a genomic bin from tumour and matched normal bam files

    Input:
        1. ref: target reference
        2. binBeg: bin begin
        3. binEnd: bin end
        4. tumourBam: indexed tumour BAM file
        5. normalBam: indexed normal BAM file
        6. confDict:
            * targetEvents   -> list with target events (INS: insertion; DEL: deletion; CLIPPING: left and right clippings; DISCORDANT: discordantly mapped read)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght

    Output:
        1. eventsDict: dictionary containing list of SV events grouped according to the SV type (only those types included in confDict[targetEvents]):
            * INS -> list of INS objects
            * DEL -> list of DEL objects
            * LEFT-CLIPPING -> list of left CLIPPING objects
            * RIGHT-CLIPPING -> list of right CLIPPING objects
            * DISCORDANT -> list of DISCORDANT objects  
    '''
    ## Search for SV events in the tumour
    eventsDict_T = collectSV(ref, binBeg, binEnd, tumourBam, confDict, 'TUMOUR')

    ## Search for SV events in the normal
    eventsDict_N = collectSV(ref, binBeg, binEnd, normalBam, confDict, 'NORMAL')

    ## Join tumour and normal lists
    eventsDict = {}

    for SV_type in eventsDict_T:        
        eventsDict[SV_type] = eventsDict_T[SV_type] + eventsDict_N[SV_type]

    return eventsDict


def collectSV(ref, binBeg, binEnd, bam, confDict, sample, supplementary = True):
    '''
    Collect structural variant (SV) events in a genomic bin from a bam file

    Input:
        1. ref: target referenge
        2. binBeg: bin begin
        3. binEnd: bin end
        4. bam: indexed BAM file
        5. confDict:
            * targetEvents   -> list with target SV events (INS: insertion; DEL: deletion; CLIPPING: left and right clippings, DISCORDANT: discordant)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght
            * overhang       -> Number of flanking base pairs around the INDEL events to be collected from the supporting read. If 'None' the complete read sequence will be collected)
        6. sample: type of sample (TUMOUR, NORMAL or None)

    Output:
        1. eventsDict: dictionary containing list of SV events grouped according to the SV type (only those types included in confDict[targetEvents]):
            * INS -> list of INS objects
            * DEL -> list of DEL objects
            * LEFT-CLIPPING -> list of left CLIPPING objects
            * RIGHT-CLIPPING -> list of right CLIPPING objects
            * DISCORDANT -> list of DISCORDANT objects  
    
    NOTE: * include secondary alignment filter???
    '''
    # Define target interval
    targetInterval = (binBeg, binEnd)

    ## Initialize dictionary to store SV events
    eventsDict = {}

    if 'INS' in confDict['targetEvents']:
        eventsDict['INS'] = []

    if 'DEL' in confDict['targetEvents']:
        eventsDict['DEL'] = []

    if 'CLIPPING' in confDict['targetEvents']:
        eventsDict['LEFT-CLIPPING'] = []
        eventsDict['RIGHT-CLIPPING'] = []
    
    if 'DISCORDANT' in confDict['targetEvents']:
        eventsDict['DISCORDANT'] = []

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, binBeg, binEnd)
    
    # For each read alignment
    for alignmentObj in iterator:

        ### 1. Filter out alignments based on different criteria:
        MAPQ = int(alignmentObj.mapping_quality) # Mapping quality

        ## Unmapped reads   
        if alignmentObj.is_unmapped == True:
            continue

        ## No query sequence available
        if alignmentObj.query_sequence == None:
            continue

        ## Aligments with MAPQ < threshold
        if (MAPQ < confDict['minMAPQ']):
            continue

        ## Duplicates filtering enabled and duplicate alignment
        if (confDict['filterDuplicates'] == True) and (alignmentObj.is_duplicate == True):
            continue
    
        # Filter supplementary alignments if FALSE. (Neccesary to avoid pick supplementary clipping reads while adding to discordant clusters in short reads mode)
        if supplementary == False and alignmentObj.is_supplementary == True:
            continue
        
        # Filtering: 
        # to do: move to a functiom
        if confDict['readFilters'] != None:
            
            # Discard alignment if mate is unmapped:
            if 'mateUnmap' in confDict['readFilters']:
                
                if alignmentObj.mate_is_unmapped:
                    continue
            
            # Filter SMS reads (reads with CIGAR #S#M#S)
            if 'SMS' in confDict['readFilters']:
                    
                firstOperation, firstOperationLen = alignmentObj.cigartuples[0]
                lastOperation, lastOperationLen = alignmentObj.cigartuples[-1]
                
                if ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation == 4) or (lastOperation == 5)):
                    continue
            
        ## 2. Collect CLIPPINGS
        if 'CLIPPING' in confDict['targetEvents']:

            left_CLIPPING, right_CLIPPING = collectCLIPPING(alignmentObj, confDict['minCLIPPINGlen'], targetInterval, sample)

            # Left CLIPPING found
            if left_CLIPPING != None:
                eventsDict['LEFT-CLIPPING'].append(left_CLIPPING)
        
            # Right CLIPPING found
            if right_CLIPPING != None:
                    
                eventsDict['RIGHT-CLIPPING'].append(right_CLIPPING)
                    
        ## 3. Collect INDELS
        if ('INS' in confDict['targetEvents']) or ('DEL' in confDict['targetEvents']):

            INDEL_events = collectINDELS(alignmentObj, confDict['targetEvents'], confDict['minINDELlen'], targetInterval, confDict['overhang'], sample)

            # Add events to the pre-existing lists                
            for INDEL_type, events in INDEL_events.items():
                eventsDict[INDEL_type] = eventsDict[INDEL_type] + events

        ## 4. Collect DISCORDANT
        if 'DISCORDANT' in confDict['targetEvents']:
            
            # Filtering: 
            # to do: move to a functiom
            if confDict['readFilters'] != None:
                
                # Discard alignment if insert size not greater than min_insertSize:
                if 'insertSize' in confDict['readFilters']:
                    
                    min_insertSize = 5000
                    insertSize = alignmentObj.template_length
                    
                    if insertSize != 0 and abs(insertSize) < min_insertSize :
                        continue
            
            DISCORDANTS = collectDISCORDANT(alignmentObj, sample)

            # Add discordant events
            for discordant in DISCORDANTS:
                eventsDict['DISCORDANT'].append(discordant)
        
    ## Close 
    bamFile.close()

    # return sv candidates
    return eventsDict


def collectCLIPPING(alignmentObj, minCLIPPINGlen, targetInterval, sample):
    '''
    For a read alignment check if the read is clipped on each side and return the corresponding clipping objects

    Input:
        1. alignmentObj: pysam read alignment object
        2. minCLIPPINGlen: minimum clipping lenght
        3. targetInterval: tuple containing begin and end position of the target genomic interval to extract events from. If 'None' all clippings will be reported
        4. sample: type of sample (TUMOUR, NORMAL or None). 

    Output:
        1. left_CLIPPING: left CLIPPING object (None if no clipping found)
        2. right_CLIPPING: right CLIPPING object (None if no clipping found)
    '''    
    # Initialize as None
    left_CLIPPING, right_CLIPPING = (None, None)

    # Determine if discordant is mate 1 or 2    
    if alignmentObj.is_read1:
        pair = '1'
    else:
        pair = '2'

    # Select first and last operation from cigar to search for clipping
    firstOperation, firstOperationLen = alignmentObj.cigartuples[0]
    lastOperation, lastOperationLen = alignmentObj.cigartuples[-1]

    ## Clipping >= X bp at the LEFT
    #  Note: soft (Operation=4) or hard clipped (Operation=5)     
    if ((firstOperation == 4) or (firstOperation == 5)) and (firstOperationLen >= minCLIPPINGlen):
        
        ## Create CLIPPING object if:
        # a) No interval specified OR 
        # b) Clipping within target interval 
        if (targetInterval == None) or (gRanges.overlap(alignmentObj.reference_start, alignmentObj.reference_start, targetInterval[0], targetInterval[1])[0]):
            
            # Create CLIPPING object
            left_CLIPPING = events.CLIPPING(alignmentObj.reference_name, alignmentObj.reference_start, alignmentObj.reference_start, firstOperationLen, 'left', pair, alignmentObj.query_name, alignmentObj.query_sequence, alignmentObj.query_alignment_start, alignmentObj, sample)

    ## Clipping > X bp at the RIGHT
    if ((lastOperation == 4) or (lastOperation == 5)) and (lastOperationLen >= minCLIPPINGlen):
 
        ## Create CLIPPING object if:
        # a) No interval specified OR 
        # b) Clipping within target interval 
        if (targetInterval == None) or (gRanges.overlap(alignmentObj.reference_end, alignmentObj.reference_end, targetInterval[0], targetInterval[1])[0]):

            # Create CLIPPING object
            right_CLIPPING = events.CLIPPING(alignmentObj.reference_name, alignmentObj.reference_end, alignmentObj.reference_end, lastOperationLen, 'right', pair, alignmentObj.query_name, alignmentObj.query_sequence, alignmentObj.query_alignment_end, alignmentObj, sample)         

    return left_CLIPPING, right_CLIPPING


def collectINDELS(alignmentObj, targetEvents, minINDELlen, targetInterval, overhang, sample):
    '''
    Collect insertions and deletions longer than a threshold that are completely spanned within an input read alignment

    Input:
        1. alignmentObj: pysam read alignment object instance
        2. targetEvents: list with target events (INS: insertion; DEL: deletion)
        3. minINDELlen: minimum INS and DEL lenght
        4. targetInterval: tuple containing begin and end position of the target genomic interval to extract events from. If 'None' all the events spanned by the read alignment will be reported
        5. overhang: Number of flanking base pairs around the SV event to be collected from the supporting read. If 'None' the complete read sequence will be collected)        
        6. sample: type of sample (TUMOUR, NORMAL or None). 

    Output:
        1. INDEL_events: dictionary containing list of SV events grouped according to the type of INDEL (only those types included in targetEvents):
            * INS -> list of INS objects
            * DEL -> list of DEL objects
    '''
    ## Initialize dict
    INDEL_events = {}

    if ('INS' in targetEvents):
        INDEL_events['INS'] = []

    if ('DEL' in targetEvents):   
        INDEL_events['DEL'] = []

    ## Initialize positions at query and ref
    posQuery = 0
    posRef = alignmentObj.reference_start

    # Iterate over the CIGAR
    for cigarTuple in alignmentObj.cigartuples:

        operation = int(cigarTuple[0])
        length = int(cigarTuple[1])
        
        ## a) INSERTION to the reference >= Xbp
        if ('INS' in targetEvents) and (operation == 1) and (length >= minINDELlen):

            ## Create INS if:
            # a) No interval specified OR 
            # b) Insertion within target interval 
            if (targetInterval == None) or (gRanges.overlap(posRef, posRef, targetInterval[0], targetInterval[1])[0]):

                # Collect piece of sequence flanking the INS event
                flankingSeq, bkpPos = (alignmentObj.query_sequence, posQuery) if overhang == None else events.pick_flanking_seq_INS(alignmentObj.query_sequence, posQuery, length, overhang)
                
                # Create INS object
                INS = events.INS(alignmentObj.reference_name, posRef, posRef, length, alignmentObj.query_name, flankingSeq, bkpPos, alignmentObj, sample)
                INDEL_events['INS'].append(INS)

        ## b) DELETION to the reference >= Xbp
        if ('DEL' in targetEvents) and (operation == 2) and (length >= minINDELlen):
            
            ## Create DEL if:
            # a) No interval specified OR 
            # b) Deletion within target interval 
            if (targetInterval == None) or (gRanges.overlap(posRef, posRef + length, targetInterval[0], targetInterval[1])[0]):

                # Collect piece of sequence flanking the DEL event
                flankingSeq, bkpPos = (alignmentObj.query_sequence, posQuery) if overhang == None else events.pick_flanking_seq_DEL(alignmentObj.query_sequence, posQuery, overhang)

                # Create DEL object
                DEL = events.DEL(alignmentObj.reference_name, posRef, posRef + length, length, alignmentObj.query_name, flankingSeq, bkpPos, alignmentObj, sample)
                INDEL_events['DEL'].append(DEL)
                
        #### Update position over reference and read sequence
        ### a) Operations consuming query and reference
        # - Op M, tag 0, alignment match (can be a sequence match or mismatch)
        # - Op =, tag 7, sequence match
        # - Op X, tag 8, sequence mismatch
        if (operation == 0) or (operation == 7) or (operation == 8):
            posQuery += length
            posRef += length

        ### b) Operations only consuming query
        # - Op I, tag 1, insertion to the reference
        # - Op S, tag 4, soft clipping (clipped sequences present in SEQ)
        elif (operation == 1) or (operation == 4):
            posQuery += length

        ### c) Operations only consuming reference
        # - Op D, tag 2, deletion from the reference
        # - Op N, tag 3, skipped region from the reference
        elif (operation == 2) or (operation == 3):
            posRef += length

        ### d) Operations not consuming query nor reference
        # - Op H, tag 5, hard clipping (clipped sequences NOT present in SEQ)
        # - Op P, tag 6, padding (silent deletion from padded reference)
        # Do not do anything

    return INDEL_events


def collectDISCORDANT(alignmentObj, sample):
    '''
    For a read alignment check if the read is discordant (not proper in pair) and return the corresponding discordant objects

    Input:
        1. alignmentObj: pysam read alignment object
        2. sample: type of sample (TUMOUR, NORMAL or None).

    Output:
        1. DISCORDANTS: list of discordant read pair events

    '''
    # Initialize discordant events list
    DISCORDANTS = []

    # If not proper pair (== discordant)
    if not alignmentObj.is_proper_pair:

        ## 1. Determine discordant orientation
        # a) Minus
        if alignmentObj.is_reverse:
            orientation = 'MINUS'

        # b) Plus
        else:
            orientation = 'PLUS'

        ## 2. Determine if discordant is mate 1 or 2
        if alignmentObj.is_read1:
            pair = '1'

        else:
            pair = '2'
                   
        ## 3. Determine number of alignment blocks
        operations = [t[0] for t in alignmentObj.cigartuples]
        nbBlocks = operations.count(3) + 1 

        ## 4. Create discordant event
        # A) Read aligning in a single block (WG or RNA-seq read no spanning a splice junction)
        if nbBlocks == 1:
            DISCORDANT = events.DISCORDANT(alignmentObj.reference_name, alignmentObj.reference_start, alignmentObj.reference_end, orientation, pair, alignmentObj.query_name, alignmentObj, sample, None)
            DISCORDANTS.append(DISCORDANT)
        
        # B) Read alignning in multiple blocks (RNA-seq read spanning one or multiple splice junctions) -> Create one discordant event per block
        else:

            blockBeg = alignmentObj.reference_start
            blockEnd = blockBeg

            # For each operation
            for cigarTuple in alignmentObj.cigartuples:

                operation = int(cigarTuple[0])
                length = int(cigarTuple[1])

                # a) End of the block -> End current block by creating a discordant event and Initiate a new block
                if operation == 3:

                    # Create discordant event for the block
                    DISCORDANT = events.DISCORDANT(alignmentObj.reference_name, blockBeg, blockEnd, orientation, pair, alignmentObj.query_name, alignmentObj, sample, None)
                    DISCORDANTS.append(DISCORDANT)

                    # Initialize new block
                    blockBeg = blockEnd + length
                    blockEnd = blockEnd + length

                # b) Extend current block
                else:
                    blockEnd = blockEnd + length   

            ## End last block by creating a discordant
            DISCORDANT = events.DISCORDANT(alignmentObj.reference_name, blockBeg, blockEnd, orientation, pair, alignmentObj.query_name, alignmentObj, sample, None)
            DISCORDANTS.append(DISCORDANT)

    return DISCORDANTS


def average_MAPQ_reads_interval(ref, beg, end, readIds, bam):
    '''
    Retrieve a set of target reads in a genomic interval. Then compute their average MAPQ

    Input:
        1. ref: target reference
        2. beg: interval begin position
        3. end: interval end position
        4. readIds: list of target read ids
        5. bam: pysam filehandler for bam file 

    Output:
        1. avMAPQ: average mapping quality
    '''
    ## 1. Collect alignments in the input interval
    iterator = bam.fetch(ref, beg, end)

    ## 2. Select only those alignments corresponding to the target reads 
    alignments = []

    for alignment in iterator:
        
        if alignment.query_name in readIds:
            alignments.append(alignment)
        
    ## 3. Compute average MAPQ for mates
    qualities = [alignment.mapping_quality for alignment in alignments]
    avMAPQ = np.mean(qualities)
    
    return avMAPQ

def samtools_index_bam(BAM, outDir):
    '''
    Index bam file using samtools

    Input:
        1. BAM: Input bam file complete path.
    
    Output:
        1. Doesn't return anything. Creates bam index files.
    '''
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    command = 'samtools index ' + BAM
    err = open(logDir + '/samtools_index_bam.err', 'w') 
    status = subprocess.call(command, stderr=err, shell=True)

    return

def BAM2FastaDict(BAM):
    '''
    Pick reads names and sequences from bam file and write tham in FASTA format.

    Input:
        1. BAM: Input bam file complete path.

    Output:
        1. fastaDict: fastaDict[readName] = readSequence
    '''

    # Read bam and store in a dictionary
    bamFile = pysam.AlignmentFile(BAM, 'rb')

    iterator = bamFile.fetch()
    
    fastaDict= {}
    # For each read alignment
    for alignmentObj in iterator:

        if alignmentObj.query_name in fastaDict.keys():
            fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
        else:
            fastaDict[alignmentObj.query_name] = []
            fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
            
    return fastaDict

def calculate_nbReads(ref, binBeg, binEnd, bam, confDict, orientation, sample, supplementary = True):
    '''
    Collect structural variant (SV) events in a genomic bin from a bam file

    Input:
        1. ref: target referenge
        2. binBeg: bin begin
        3. binEnd: bin end
        4. bam: indexed BAM file
        5. confDict:
            * targetEvents       -> list with target SV events (INS: insertion; DEL: deletion; CLIPPING: left and right clippings, DISCORDANT: discordant)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght
            * overhang       -> Number of flanking base pairs around the INDEL events to be collected from the supporting read. If 'None' the complete read sequence will be collected)
        6. sample: type of sample (TUMOUR, NORMAL or None)
    Output:
        1. NB_READ: number of alignments in interval. All included
    '''
    # Define target interval
    targetInterval = (binBeg, binEnd)
    
    ## Initialize counter
    NB_READ = 0
        
    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")
    
    ## Extract alignments
    iterator = bamFile.fetch(ref, binBeg, binEnd)
    
    # For each read alignment
    for alignmentObj in iterator:
        
        ### 1. Filter out alignments based on different criteria:
        MAPQ = int(alignmentObj.mapping_quality) # Mapping quality
        
        ## Unmapped reads   
        if alignmentObj.is_unmapped == True:
            continue
        
        ## No query sequence available
        if alignmentObj.query_sequence == None:
            continue
        
        ## Aligments with MAPQ < threshold
        if (MAPQ < confDict['minMAPQ']):
            continue
        
        ## Duplicates filtering enabled and duplicate alignment
        if (confDict['filterDuplicates'] == True) and (alignmentObj.is_duplicate == True):
            continue
        
        # Filter supplementary alignments if FALSE. (Neccesary to avoid pick supplementary clipping reads while adding to discordant clusters in short reads mode)
        if supplementary == False and alignmentObj.is_supplementary == True:
            continue
        
        # Filtering: 
        # to do: move to a functiom
        if confDict['readFilters'] != None:
            
            # Discard alignment if mate is unmapped:
            if 'mateUnmap' in confDict['readFilters']:
                
                if alignmentObj.mate_is_unmapped:
                    continue
            
            # Filter SMS reads (reads with CIGAR #S#M#S)
            if 'SMS' in confDict['readFilters']:
                    
                firstOperation, firstOperationLen = alignmentObj.cigartuples[0]
                lastOperation, lastOperationLen = alignmentObj.cigartuples[-1]
                
                if ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation == 4) or (lastOperation == 5)):
                    continue
        
        if orientation == "+" and alignmentObj.is_reverse:
            continue
        
        if orientation == "-" and not alignmentObj.is_reverse:
            continue
             
        NB_READ += 1
        
    ## Close 
    bamFile.close()
    
    # return sv candidates
    return NB_READ