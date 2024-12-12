'''
Module 'annotation' - Contains functions for the annotation of genomic intervals according to different annotation resources
'''

## DEPENDENCIES ##
# External
import os
import subprocess
from operator import itemgetter
import pickle

# Internal
from GAPI import unix
from GAPI import formats
from GAPI import databases
from GAPI import log
from GAPI import gRanges


def newestFile(path):
    '''
    Get the most recently modified file in path
    
    Input:
    1. path: directory to be interrogated
    
    Output:
    2. newestPath: newest path or file in the indicated directory
    '''
    files = os.listdir(path)
    paths = [os.path.join(path, basename) for basename in files]
    newestPath = max(paths, key=os.path.getctime)
    
    return newestPath

def load_annotations(annotations2load, refLengths, annotationsDir, germlineMEI, threads, outDir, tdEnds = ['3']):
    '''
    Load a set of annotation files in bed formats into a bin database

    Input:
        1. annotations2load: list of annotations to load. Annotations available: REPEATS, TRANSDUCTIONS and EXONS
        2. refLengths: Dictionary containing reference ids as keys and as values the length for each reference  
        3. annotationsDir: Directory containing annotation files
        4. germlineMEI: Bed file containing set of known germline MEI. None if not available
        5. threads: number of threads used to parallelize the bin database creation
        6. outDir: Output directory
        7. tdEnds: Default [3]. [3, 5]
    
    Output:
        1. annotations: dictionary containing one key per type of annotation loaded and bin databases containing annotated features as values (None for those annotations not loaded)
    '''
    
    ## 1. Check for cached dictionary
    cachePath = annotationsDir + '/.cache'
    cacheFile = cachePath + '/' + '_'.join(annotations2load) + '.pkl'
    
    # If cache exist and it's the most recent file in dir
    if os.path.exists(cacheFile) and newestFile(annotationsDir) == cachePath: 
            
            # Load annotations from cache
            pickle_in = open(cacheFile, "rb")
            annotations = pickle.load(pickle_in)
            
            return annotations
    
    ## 2. Initialize dictionary
    annotations = {}
    annotations['REPEATS'] = None
    annotations['REPEATS-POLYA'] = None
    annotations['RETROTRANSPOSONS'] = None
    annotations['TRANSDUCTIONS'] = None
    annotations['EXONS'] = None
    annotations['GERMLINE-MEI'] = None

    ## Create output directory
    unix.mkdir(outDir)

    ## 3. Load annotated repeats into a bin database
    if 'REPEATS' in annotations2load:

        repeatsBed = annotationsDir + '/repeats.bed'
        annotations['REPEATS'] = formats.bed2binDb(repeatsBed, refLengths, threads)
        
    elif 'REPEATS-L1' in annotations2load:
    
        repeatsBed = annotationsDir + '/repeats.L1.pA.bed'
        annotations['REPEATS'] = formats.bed2binDb(repeatsBed, refLengths, threads)

    ## 4. Load annotated repeats into a bin database
    if 'RETROTRANSPOSONS' in annotations2load:

        repeatsBed = annotationsDir + '/retrotransposons.bed'
        annotations['RETROTRANSPOSONS'] = formats.bed2binDb(repeatsBed, refLengths, threads)
        
    ## 5. Create germline MEI database
    if 'REPEATS-POLYA' in annotations2load:
        
        polyABed = annotationsDir + '/polyA.bed'
        annotations['REPEATS-POLYA'] = formats.bed2binDb(polyABed, refLengths, threads)
        
    ## 6. Create transduced regions database
    if 'TRANSDUCTIONS' in annotations2load:

        ## Create bed file containing transduced regions
        sourceBed = annotationsDir + '/srcElements.bed'
        # buffer equals -150 to avoid the end of the src element
        transducedPath = databases.create_transduced_bed(sourceBed, tdEnds, 10000, -150, outDir)
        
        ## Load transduced regions into a bin database
        annotations['TRANSDUCTIONS'] = formats.bed2binDb(transducedPath, refLengths, threads)

    ## 7. Create exons database
    if 'EXONS' in annotations2load:

        exonsBed = annotationsDir + '/exons.bed'
        annotations['EXONS'] = formats.bed2binDb(exonsBed, refLengths, threads)

    ## 8. Create germline MEI database
    if 'GERMLINE-MEI' in annotations2load and germlineMEI:
        
        annotations['GERMLINE-MEI'] = formats.bed2binDb(germlineMEI, refLengths, threads)

    ## 9. Create the cache if does not exist
    unix.mkdir(cachePath)
    pickleOut = open(cacheFile, "wb")
    pickle.dump(annotations, pickleOut)
    pickleOut.close()
        
    return annotations

def annotate(events, steps, annotations, annovarDir, outDir):
    '''
    Annotate each input event interval based on different annotation resources.

    Input: 
        1. events: List containing input events to be annotated. Events should be objects containing ref, beg and end attributes.
        2. steps: List containing annotation steps to be performed. Possible values: REPEAT and GENE.
        3. annotations: Dictionary containing bin databases for annotations.  
        4. annovarDir: Directory containing annovar reference databases. 'None' if gene annotation no enabled 
        5. outDir: Output directory

    Output:
        1. New 'repeatAnnot' attribute set for each input event. 
        'repeatAnnot' is a list of dictionaries. Each dictionary contains information pertaining to one overlapping repeat

        2. New 'geneAnnot' attribute set for each input event. 
        'geneAnnot' is a tuple(region,gene) 
    '''
    ## 1. Perform repeat based annotation if enabled
    msg = '1. Repeat based annotation if enabled'
    log.subHeader(msg)   

    if 'REPEAT' in steps:

        msg = 'Perform repeats annotation'
        log.info(msg)   
        repeats_annotation(events, annotations['REPEATS'], 200)

    ## 2. Perform gene-based annotation if enabled
    msg = '2. Gene-based annotation if enabled'
    log.subHeader(msg)

    if 'GENE' in steps:
        msg = 'Perform gene-based annotation'
        log.info(msg)  
        gene_annotation(events, annovarDir, outDir)

    ## 3. Perform known germline MEI annotation if enabled
    msg = '3. Perform known germline MEI annotation if enabled'
    log.subHeader(msg)

    if 'GERMLINE-MEI' in steps:
        msg = 'Perform known germline MEI annotation if enabled'
        log.info(msg)  
        germline_MEI_annotation(events, annotations['GERMLINE-MEI'], 150)


def annotate_interval(ref, beg, end, annotDb):
    '''
    Intersect input interval (ref:beg-end) with a given annotation 

    Input: 
        1. ref: reference id
        2. beg: begin position
        3. end: end position
        4. annotDb: dictionary containing annotated features organized per chromosome (keys) into genomic bins (values)

    Output:
        1. sortedOverlaps. List of lists sorted in decreasing percentage of overlap. Each tuple corresponds to one overlapping event and is composed by 3 elements: 
            1. Overlapping event
            2. Number of overlapping base pairs
            3. Percentage of base pairs of the input interval that are overlapping  
            4. Tuple with input interval coordinates overlapping with the event
    '''

    # a) Annotated features available in the same ref 
    if ref in annotDb:
            
        ## Select features bin database for the corresponding reference 
        binDb = annotDb[ref]        

        ## Retrieve all the annotated features overlapping with the input interval
        overlaps = binDb.collect_interval(beg, end, 'ALL')    

        ## Order overlapping features in decreasing order of perc of overlap
        sortedOverlaps = sorted(overlaps, key=lambda x: x[2], reverse=True)
         
    # b) No feature in the same ref as the interval
    else:
        sortedOverlaps = []

    return sortedOverlaps
        
def repeats_annotation(events, repeatsDb, buffer):
    '''
    For each input event assess if overlaps with an annotated repeat in the reference genome

    Input: 
        1. events: list containing input events to be annotated. Events should be objects containing ref, beg and end attributes.
        2. repeatsDb: dictionary containing annotated repeats organized per chromosome (keys) into genomic bins (values)
        3. buffer: number of base pairs to extend begin and end coordinates for each event prior assessing overlap

    Output:
        New 'repeatAnnot' attribute set for each input event. 
        'repeatAnnot' is a list of dictionaries. Each dictionary contains information pertaining to one overlapping repeat
    '''

    ## Assess for each input event if it overlaps with an annotated repeat
    for event in events:

        # A) Annotated repeat in the same ref where the event is located
        if event.ref in repeatsDb:
            
            ### Select repeats bin database for the corresponding reference 
            repeatsBinDb = repeatsDb[event.ref]        

            ### Retrieve all the annotated repeats overlapping with the event interval
            overlaps = repeatsBinDb.collect_interval(event.beg - buffer, event.end + buffer, 'ALL')    

            ### Compute distance between the annotated repeat and the raw interval
            annotatedRepeats = []

            ## For each intersection
            for overlap in overlaps:

                repeat = {}
                repeat['family'] = overlap[0].optional['family']
                repeat['subfamily'] = overlap[0].optional['subfamily']

                overlapLen = overlap[1] 
                boolean = gRanges.overlap(event.beg, event.end, overlap[0].beg, overlap[0].end)[0]

                # a) Overlapping raw intervals, distance == 0 
                if boolean:
                    distance = 0

                # b) Not overlapping raw intervals. Compute distance
                else:
                    distance = abs(overlapLen - buffer)

                repeat['distance'] = distance
                annotatedRepeats.append(repeat)

        # B) No repeat in the same ref as the event
        else:
            annotatedRepeats = []
        
        ## Add repeat annotation as attribute 
        event.repeatAnnot = annotatedRepeats


def gene_annotation(events, annovarDir, outDir):
    '''
    Perform gene-based annotation for a list of input events
 
    Input: 
        1. events: List containing input events to be annotated. Events should be objects containing ref, beg and end attributes.
        2. annovarDir: Directory containing the two files used by ANNOVAR to perform gene based-annotation:
                a) build_annot.txt     - Text file containing annotated transcript coordinates
                b) build_annotMrna.fa  - Fasta containing annotated transcript sequences
        3. outDir: Output directory

    Output:
    
        New 'geneAnnot' attribute set for each input event. 
        'geneAnnot' is a tuple(region,gene) 
    '''

    ## 1. Create output directory
    unix.mkdir(outDir)

    ## 2. Create input file containing events intervals for ANNOVAR 
    create_annovar_input(events, 'events.annovar', outDir)
    annovarInput = outDir + '/events.annovar'

    ## 3. Annotate events intervals with ANNOVAR 
    out1, out2 = run_annovar(annovarInput, annovarDir, outDir)

    ## 4. Add gene annotation info to the events
    addGnAnnot2events(events, out1)
    
    ## Do cleanup
    unix.rm([annovarInput, out1, out2])

def addGnAnnot2events(events, out1):
    '''
    Read annovar output file and incorporate gene annotation information to the corresponding event objects

    Input: 
        1. events: List containing input events to be annotated. Events should be objects containing ref, beg and end attributes.
        2. out1: Annovar output file 1 (region annotation for all the variants) 

    Output:
        New 'geneAnnot' attribute set for each input event. 
        'geneAnnot' is a tuple(region,gene) 
    '''
    ## 1. Organize events into a dict
    eventsDict = {}

    for event in events:
        
        beg = event.refLeftBkp if event.refLeftBkp is not None else event.beg
        end = event.refRightBkp if event.refRightBkp is not None else event.end
        
        name = event.ref + ':' + str(beg) + '-' + str(end)
        eventsDict[name] = event

    ## 2. Add to each event gene annotation info    
    out1File = open(out1, "r")

    # Read line by line adding the relevant info to the corresponding event in each iteration
    for line in out1File:

        fields = line.split()
        region = fields[0]
        gene = fields[1]
        name = fields[8]             
        eventsDict[name].geneAnnot = (region, gene)
        
def germline_MEI_annotation(events, MEIDb, buffer):
    '''
    For each input event assess if overlaps with an already known germline MEI polymorphism

    Input: 
        1. events: list containing input events to be annotated. Events should be objects containing ref, beg and end attributes.
        2. MEIDb: dictionary containing known germline MEI organized per chromosome (keys) into genomic bin databases (values)
        3. buffer: number of base pairs to extend begin and end coordinates for each event prior assessing overlap

    Output:
        New 'germlineDb' attribute set for each input event. Attribute is a string
        with the germline MEI databases where the MEI has already been reported
    '''
    ## Assess for each input event if it overlaps with an germline MEI repeat
    for event in events:

        # Skip event if family not available or no MEI on the same reference
        if ('FAMILY' not in event.SV_features) or (event.ref not in MEIDb):
            event.germlineDb = None            
            continue

        ### Select bin database for the corresponding reference 
        binDb = MEIDb[event.ref]        

        ### Retrieve overlapping known germline MEI if any
        overlaps = binDb.collect_interval(event.beg - buffer, event.end + buffer, event.SV_features['FAMILY']) 

        ### a) Known germline MEI overlapping the event
        if overlaps:
            event.germlineDb = overlaps[0][0].optional['database']                     

        ### b) No overlap found                    
        else:
            event.germlineDb = None

def create_annovar_input(events, fileName, outDir):
    '''
    Write events intervals into a format compatible with annovar 

    Input:
        1. events: List containing input events. Events should be objects containing ref, beg and end attributes.
        2. fileName: Output file name 
        3. outDir: Output directory
    '''
    ## 1. Write header
    outPath = outDir + '/' + fileName
    outFile = open(outPath, 'w')

    ## 2. Write events
    # For each event
    for event in events:
        
        beg = event.refLeftBkp if event.refLeftBkp is not None else event.beg
        end = event.refRightBkp if event.refRightBkp is not None else event.end

        name = event.ref + ':' + str(beg) + '-' + str(end)
        fields = [event.ref, str(beg), str(end), '0', '0', 'comments: ' + name]
        row = "\t".join(fields)
        
        # Add entry to BED
        outFile.write(row + '\n')


def read_annovar_out1(out1):
    '''
    Read annovar out1 and store info in a dictionary

    Input:
        1. out1: Annovar output file 1 (region annotation for all the variants) 

    Output:
        1. out1Dict: Dictionary with comment field as keys and tuples(region,gene) as values        
    '''

    ## Open annovar out file
    with open(out1, "r") as out1File:
    
        out1Dict = {}

        # Read line by line adding the relevant info to the dict in each iteration
        for line in out1File:
            fields = line.split()
            region = fields[0]
            gene = fields[1]
            name = fields[8] 
            
            out1Dict[name] = (region, gene)

    return out1Dict
            

def run_annovar(inputFile, annovarDir, outDir):
    '''
    For each input entry perform gene-based annotation with Annovar

    Input: 
        1. inputFile: Annovar-like file containing regions to be annotated
        2. annovarDir: Directory where annovar annotation files are located
        3. outDir: Output directory

    Output: 
        1. out1: Annovar output file 1 (region annotation for all the variants) 
        2. out2: Annovar output file 2 (amino acid changes as a result of the exonic variant) 

    NOTE: perl annotate_variation.pl call should be available as 'ANNOVAR' environmental variable
    NOTE: check annovar documentation for explanation about how to interpret output files:

    http://annovar.openbioinformatics.org/en/latest/user-guide/gene/#output-file-1-refseq-gene-annotation
    http://annovar.openbioinformatics.org/en/latest/user-guide/gene/#output-file-2-refseq-gene-annotation
    '''

    ## Run annovar in the input file
    ANNOVAR = os.environ['ANNOVAR']
    out = open(outDir + '/annovar.out', 'w') 
    err = open(outDir + '/annovar.err', 'w') 
    command = ANNOVAR + ' -buildver build -out ' + outDir + '/annovar -dbtype annot ' + inputFile + ' ' + annovarDir
    status = subprocess.call(command, stdout=out, stderr=err, shell=True)

    ## Return results
    out1 = outDir + '/annovar.variant_function'
    out2 = outDir + '/annovar.exonic_variant_function'

    return out1, out2


def intersect_mate_annotation(discordants, annotation, targetFields, readSize):
    '''
    For each input read assess if the mate aligns over an annotated feature

    Input: 
        1) discordants: list containing input discordant read pair events
        2) annotation: dictionary containing annotated features organized per chromosome (keys) into genomic bins (values)
        3) targetFields: List of optional fields to be used to determine overlapping feature name
        4) readSize: read size

    Output:
        1) matesIdentity: dictionary containing lists of discordant read pairs organized taking into account their orientation and if the mate aligns in an annotated feature 
                               This info is encoded in the dictionary keys as follows. Keys composed by 3 elements separated by '-':
                                
                                    - Orientation: read orientation (PLUS or MINUS)
                                    - Event type: DISCORDANT   
                                    - featureType: Feature type or 'None' if mate does not align in a retrotransposon
    '''
    matesIdentity = {}

    ##  For each input discordant intersect mate alignment coordinates with the provided annotation 
    for discordant in discordants:
        
        # A) Annotated feature in the same ref where the mate aligns
        if discordant.mateRef in annotation:

            ## Select features bin database for the corresponding reference 
            featureBinDb = annotation[discordant.mateRef]        
 
            ## Retrieve all the annotated features overlapping with the mate alignment interval
            buffer = readSize if readSize else 100
            overlappingFeatures = featureBinDb.collect_interval(discordant.mateStart, discordant.mateStart + buffer, 'ALL')
            
            ## Determine mate status 
            # a) Mate does not align within an annotated feature
            if len(overlappingFeatures) == 0:
                featureType = 'None'

            # b) Mate aligns within a single feature
            elif len(overlappingFeatures) == 1:
                featureType = '_'.join([overlappingFeatures[0][0].optional[targetField] for targetField in targetFields])

            # c) Mate overlaps multiple features
            else:
                overlappingFeatures = sorted(overlappingFeatures, key=itemgetter(1), reverse=True)
                featureType = '_'.join([overlappingFeatures[0][0].optional[targetField] for targetField in targetFields])
                
        # B) No feature in the same ref as mate
        else:
            featureType = 'None'

        ## Set discordant identity
        discordant.identity = featureType

        ## Add discordant read pair to the dictionary
        identity = discordant.orientation + '_DISCORDANT_' + featureType

        # a) There are already discordant read pairs with this identity
        if identity in matesIdentity:
            matesIdentity[identity].append(discordant)

        # b) First discordant read pair with this identity
        else:
            matesIdentity[identity] = [ discordant ] 
    
    return matesIdentity
