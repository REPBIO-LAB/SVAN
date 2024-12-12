'''
Module 'clusters' - Contains classes for dealing with genomic variation
'''

## DEPENDENCIES ##
# External
import sys
import multiprocessing as mp
import numpy as np
import pysam
import math
import itertools
import os
from operator import itemgetter
from collections import Counter
import scipy
import time

# Internal
from GAPI import log
from GAPI import formats
from GAPI import unix
from GAPI import clustering
from GAPI import events
from GAPI import structures
from GAPI import alignment
from GAPI import bamtools 
from GAPI import assembly
from GAPI import repeats
from GAPI import retrotransposons
from GAPI import filters
from GAPI import bkp
from GAPI import annotation
from GAPI import sequences
from GAPI import gRanges
from GAPI import virus
from GAPI import stats

###############
## FUNCTIONS ##
###############

def create_cluster(events, clusterType):
    '''
    Function to create a cluster object instance

    Input:
        1. events: list of events/clusters that will compose the cluster
        2. clusterType: type of cluster (META: metacluster; INS: insertion; DEL: deletion; CLIPPING: clipping; DISCORDANT: discordant paired-end)

    Output:
        1. cluster: cluster object instance
    '''
    cluster = ''

    ## a) Create META cluster
    # Note: events should correspond to a list of clusters 
    if (clusterType == 'META'):
        cluster = META_cluster(events)

    ## b) Create INS cluster
    elif (clusterType == 'INS'):
        cluster = INS_cluster(events)

    ## c) Create DEL cluster
    elif (clusterType == 'DEL'):
        cluster = DEL_cluster(events)

    ## d) Create CLIPPING cluster
    elif (clusterType == 'CLIPPING') or (clusterType == 'LEFT-CLIPPING') or (clusterType == 'RIGHT-CLIPPING'):
        cluster = CLIPPING_cluster(events)

    ## e) Create DISCORDANT cluster
    elif 'DISCORDANT' in clusterType:
        cluster = DISCORDANT_cluster(events)

    ## f) Create SUPPLEMENTARY cluster 
    elif 'SUPPLEMENTARY' in clusterType:
        cluster = SUPPLEMENTARY_cluster(events)

    ## g) Create INS_VCF cluster
    elif 'INS_VCF' in clusterType:
        cluster = formats.INS_cluster(events)

    ## h) Unexpected cluster type
    else:
        log.info('Error at \'create_cluster\'. Unexpected cluster type')
        sys.exit(1)

    return cluster


def merge_clusters(clusters, clusterType):
    '''
    Merge a set of clusters/metaclusters into a single cluster/metacluster instance

    Input:
        1. clusters: list of clusters/metaclusters that will be merged
        2. clusterType: type of cluster to be merged (INS: insertion; DEL: deletion; CLIPPING: clipping; META: metacluster)

    Output:
        1. cluster: merged cluster/metacluster instance
    '''
    # A) Merge metaclusters
    if clusterType == 'META':
        subclusters = []
         
        for metacluster in clusters:
            # NOTE MERGE SR2020: To avoid key META error in clustering.reciprocal_overlap_clustering
            #subclusters = subclusters + list(metacluster.subclusters.values())
            subclusters = subclusters + list(metacluster.rawSubclusters)

        mergedCluster = create_cluster(subclusters, clusterType)

        # NOTE MERGE SR2020: To avoid key META error in clustering.reciprocal_overlap_clustering
        for metacluster in clusters:
            for clusterNew in metacluster.rawSubclusters:
                clusterNew.clusterId = mergedCluster.id

    # B) Merge standard clusters
    else:

        events = []
        for cluster in clusters:
            events = events + cluster.events

        mergedCluster = create_cluster(events, clusterType)

    return mergedCluster


def cluster_by_matePos(discordants, refLengths, minClusterSize):
    '''
    Apply an extra clustering step to discordant read pair clusters based on mate position

    Input:
        1. discordants: list of discordant clusters
        2. refLengths: dictionary containing references as keys and their lengths as values
        3. minClusterSize: minimum cluster size

    Output:
        1. outDiscordants: list of discordant clusters after applying mate position based clustering
    '''     
    outDiscordants = []

    ## For each cluster
    for cluster in discordants:

        ## Cluster by mate position
        newClusters = cluster_events_by_matePos(cluster.events, refLengths, minClusterSize)
            
        ## Add newly created clusters to the list
        outDiscordants = outDiscordants + newClusters

    return outDiscordants


def cluster_by_identity(clusters, refLengths, minClusterSize):
    '''
    Apply an extra clustering step to clusters based on mate or SA identity

    Input:
        1. clusters: list of clusters
        2. refLengths: dictionary containing references as keys and their lengths as values
        3. minClusterSize: minimum cluster size

    Output:
        1. outClusters: list of clusters after applying identity based clustering
    '''     
    outClusters = []

    ## For each cluster
    for cluster in clusters:

        ## Cluster by mate position
        newClusters = cluster_events_by_identity(cluster.events, refLengths, minClusterSize)
            
        ## Add newly created clusters to the list
        outClusters = outClusters + newClusters

    return outClusters


def cluster_by_supplPos(clippings, refLengths, minClusterSize, clippingSide):
    '''
    Apply an extra clustering step to clipping clusters based on suppl. alignment position

    Input:
        1. clippings: list of clipping clusters
        2. refLengths: dictionary containing references as keys and their lengths as values
        3. minClusterSize: minimum cluster size
        4. clippingSide: 'LEFT-CLIPPING' or 'RIGHT-CLIPPING'

    Output:
        1. outClippings: list of clipping clusters after applying supplementary alignment position based clustering
    '''     
    outClippings = []

    ## For each cluster
    for cluster in clippings:

        ## Perform extra clustering based on suppl. alignment positions
        newClusters = cluster_events_by_supplPos(cluster, refLengths, minClusterSize, clippingSide)
            
        ## Add newly created clusters to the list
        outClippings = outClippings + newClusters

    return outClippings

def cluster_events_by_supplPos(clippingCluster, refLengths, minClusterSize, clippingSide):
    '''
    Cluster events in an input clipping cluster according to the position of their supplementary alignments
    
    Input:
        1. clippingCluster: clipping cluster
        2. refLengths: dictionary containing references as keys and their lengths as values
        3. minClusterSize: minimum cluster size
        4. clippingSide: 'LEFT-CLIPPING' or 'RIGHT-CLIPPING'

    Output:
        1. clippingClusters: list of clipping clusters
    '''   
    ## 1. Organize clippings at the cluster into a dictionary 
    clippingsDict = {}

    for clipping in clippingCluster.events:
        clippingsDict[clipping.id] = clipping
    
    ## 2. Retrieve complementary suppl. alignments for each clipping event 
    complementaryAlignments = clippingCluster.search4complAlignments()

    ## 3. Cluster supplementary alignment positions
    supplClusters = clippingCluster.cluster_suppl_positions(complementaryAlignments)
    
    ## 4. Create new clipping clusters of discordants based on mate clusters
    clippingClusters = []

    supplClusters = structures.dict2list(supplClusters)
 
    # For each mate cluster in the reference
    for supplCluster in supplClusters:
        
        # Retrieve original discordants for mates
        clippings = [clippingsDict[supplAlign.clippingId] for supplAlign in supplCluster.events]

        # Create cluster
        clippingCluster = create_cluster(clippings, clippingSide)

        # Add supplementary alignments cluster to the clipping cluster
        clippingCluster.supplCluster = supplCluster

        # Add clipping cluster to the list
        clippingClusters.append(clippingCluster)

    return clippingClusters

def cluster_events_by_matePos(discordants, refLengths, minClusterSize):
    '''
    Cluster discordant read pair events based on their mate alignment position

    Input:
        1. discordants: list of discordant events
        2. refLengths: dictionary containing references as keys and their lengths as values
        3. minClusterSize: minimum cluster size

    Output:
        1. discordantClusters: list of discordant clusters
    '''   
    ## 1. Organize discordant into a dictionary according to supporting read id
    discordantsDict = {}

    for discordant in discordants:
        discordantsDict[discordant.fullReadName()] = discordant
    
    ## 2. Produce discordant objects for mates:
    mates = events.discordants2mates(discordants)

    ## 3. Organize mates into a bin database prior clustering
    ## Create dictionary containing mates
    matesDict = events.events2nestedDict(mates, 'DISCORDANT_MATE')

    ## Create bin database 
    matesBinDb = structures.create_bin_database(refLengths, matesDict)

    ## 4. Cluster mates according to their alignment positions
    mateClusters = []

    # For each reference
    for ref in matesBinDb:
        binDb = matesBinDb[ref]
        binLevel = binDb.binSizes[0]
        # clusters = clustering.distance_clustering(binDb, binLevel, ['DISCORDANT_MATE'], 'DISCORDANT', binLevel, minClusterSize)
        clusters = clustering.distance_clustering(binDb, binLevel, ['DISCORDANT_MATE'], 'DISCORDANT', 500, minClusterSize)
        mateClusters = mateClusters + clusters
    
    ## 5. Make clusters of discordants based on mate clusters
    discordantClusters = []

    # For each mate cluster in the reference
    for mateCluster in mateClusters:
        
        # Retrieve original discordants for mates
        discordants = [discordantsDict[mate.fullReadName_mate()] for mate in mateCluster.events]

        # Create cluster
        discordantCluster = create_cluster(discordants, 'DISCORDANT')

        # Add group to the list
        discordantClusters.append(discordantCluster)

    return discordantClusters

def create_clusters(eventsBinDb, confDict):
    '''
    Group SV events into distinct clusters according to their SV type
    
    Input:
        1. eventsBinDb: Data structure containing a set of events organized in genomic bins
        2. confDict: 
            * maxInsDist: Maximum distance bewteen two adjacent INS to be clustered together
            * maxBkpDist: Maximum distance bewteen two adjacent breakpoints for CLIPPING clustering    
            * minClusterSize: minimum number of reads composing a root cluster
            * minPercOverlap: minimum percentage of reciprocal overlap for DEL clustering 

    Output:
        1. clustersBinDb: bin database structure containing SV clusters
    ''' 

    ## 1. Create clusters ##
    clustersDict = {}

    ## For each event type, group events into clusters 
    for SV_type in eventsBinDb.eventTypes:

        ## A) Perfom clustering for INS events
        if SV_type == 'INS':

            binLevel = eventsBinDb.binSizes[0]
            clustersDict[SV_type] = clustering.distance_clustering(eventsBinDb, binLevel, [SV_type], SV_type, confDict['maxInsDist'], confDict['minClusterSize'])      

        ## B) Perform clustering for CLIPPING events
        elif SV_type in ['CLIPPING', 'RIGHT-CLIPPING', 'LEFT-CLIPPING']:

            binLevel = eventsBinDb.binSizes[0]
            clustersDict[SV_type] = clustering.distance_clustering(eventsBinDb, binLevel, [SV_type], SV_type, confDict['maxBkpDist'], confDict['minClusterSize'])      

        ## C) Perform clustering based on reciprocal overlap for DEL 
        elif SV_type == 'DEL':

            clustersDict[SV_type] = clustering.reciprocal_overlap_clustering(eventsBinDb, 20, confDict['minClusterSize'], [SV_type], 0, SV_type)

        ## D) Perform clustering based on reciprocal overlap for DISCORDANT
        elif SV_type == 'DISCORDANT':

            clustersDict[SV_type] = clustering.reciprocal_overlap_clustering(eventsBinDb, 1, confDict['minClusterSize'], [SV_type], confDict['equalOrientBuffer'], SV_type)    

    ## 2. Organize clusters into bins ##    
    binSizes = [100, 1000, 10000, 100000, 1000000]
    clustersBinDb = structures.create_bin_database_interval(eventsBinDb.ref, eventsBinDb.beg, eventsBinDb.end, clustersDict, binSizes)

    return clustersBinDb


def polish_clusters(clustersBinDb, minClusterSize):
    '''
    Apply set of steps for refining clusters of SV events

    Input:
        1. clustersBinDb: Data structure containing a set of clusters organized in genomic bins  
        2. minClusterSize: minimum cluster size

    Output: 
        It does not return variables as output. Just modify the bin database
    '''
    ## 1. Polish INS clusters
    if 'INS' in clustersBinDb.eventTypes:
        
        # Collect all clusters
        INS_clusters = clustersBinDb.collect(['INS'])

        # Correct INS fragmentation
        for cluster in INS_clusters:    

            cluster.correct_fragmentation()

        # Search for subclusters and remove outliers
        for cluster in INS_clusters:

            subclusters = cluster.identify_subclusters(minClusterSize)

            # Replace cluster by newly created subclusters in the bin database
            if subclusters:

                clustersBinDb.remove([cluster], 'INS')
                clustersBinDb.add(subclusters, 'INS')

    ## 2. Polish DEL clusters (TO DO) 

    ## 3. Polish CLIPPING clusters (TO DO)


def create_metaclusters(clustersBinDb, buffer = 200):
    '''    
    Group SV cluster events into metaclusters

    Input:
        1. clustersBinDb: Data structure containing a set of clusters organized in genomic bins
        2. buffer: buffer to extend cluster coordinates

    Output:
        1. metaclusters: list containing newly created metaclusters
    '''
    metaclusters = clustering.reciprocal_overlap_clustering(clustersBinDb, 1, 1, clustersBinDb.eventTypes, buffer, 'META')

    return metaclusters


def create_metaclusters_distanceClustering(clustersBinDb, binSize, buffer):
    '''    
    Group SV cluster events into metaclusters

    Input:
        1. clustersBinDb: Data structure containing a set of clusters organized in genomic bins
        2. buffer: buffer to extend cluster coordinates

    Output:
        1. metaclusters: list containing newly created metaclusters
    '''
    metaclusters = clustering.distance_clustering(clustersBinDb, binSize, clustersBinDb.eventTypes, 'META', buffer, 1)
    
    return metaclusters

def SV_type_metaclusters(metaclusters, minINDELlen, technology, rootOutDir):
    '''
    Infer the SV type supported by each metacluster

    Input:
        1. metaclusters: list of metaclusters
        2. minINDELlen: minimum INS and DEL lenght
        3. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        4. rootOutDir: root output directory

    Output:
        1. metaclustersSVType: dictionary containing one key per SV type and the list of metaclusters identified as value
    '''
    metaclustersSVType = {}

    for metacluster in metaclusters:

        metaInterval = '_'.join([str(metacluster.ref), str(metacluster.beg), str(metacluster.end)])
        outDir = rootOutDir + '/' + metaInterval
        metacluster.determine_SV_type(minINDELlen, technology, outDir)

        # A) Initialize list containing metaclusters of a given SV type
        if metacluster.SV_type not in metaclustersSVType:
            metaclustersSVType[metacluster.SV_type] = [metacluster]

        # B) Add metacluster to the list        
        else:
            metaclustersSVType[metacluster.SV_type].append(metacluster)    

    return metaclustersSVType


def create_consensus(metaclusters, confDict, reference, targetSV, rootOutDir):
    '''
    Generate consensus events for a set of metacluster objects

    Input:
        1. metaclusters: Dictionary containing one key per SV type and the list of metaclusters identified as value
        2. confDict: 
            * technology     -> sequencing technology (NANOPORE, PACBIO or ILLUMINA)
            * rounds         -> number of polishing rounds to be attempled. 0 means no polishing
            * targetSV       -> list with target SV (INS: insertion; DEL: deletion; CLIPPING: left and right clippings)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght
            * overhang       -> Number of flanking base pairs around the INDEL events to be collected from the supporting read. If 'None' the complete read sequence will be collected)

        3. reference: Path to reference genome in fasta format    
        4. targetSV: Target SV types to generate consensus
        5. rootOutDir: Root output directory
    ''' 

    ## For each type of SV 
    for SV in targetSV:

        ## Abort if no metacluster from this SV type has been identified
        if SV not in metaclusters:
            continue

        ## For each metacluster
        for metacluster in metaclusters[SV]:
            metaInterval = '_'.join([str(metacluster.ref), str(metacluster.beg), str(metacluster.end)])
            outDir = rootOutDir + '/' + metaInterval

            ## 1. Polish metacluster´s consensus sequence
            metacluster.polish(confDict, reference, outDir)

            ## 2. Obtain consensus metacluster´s event
            metacluster.consensus_event(confDict, reference, 10000, outDir)

            ## Cleanup
            unix.rm([outDir])
    

def lighten_up_metaclusters(metaclusters):
    '''
    Make metacluster objects lighter by removing events and subcluster objects.

    Collect relevant info for downstream analysis as metacluster attributes

    Input:
        1. metaclusters: Dictionary containing one key per SV type and the list of metaclusters identified as value
    '''
    # For each SV type
    for SV_type in metaclusters:

        # For each metacluster
        for metacluster in metaclusters[SV_type]:

            ## Set some object attributes before lightening up
            metacluster.nbTotal, metacluster.nbTumour, metacluster.nbNormal, metacluster.nbINS, metacluster.nbDEL, metacluster.nbCLIPPING = metacluster.nbEvents()
            metacluster.nbReadsTotal, metacluster.nbReadsTumour, metacluster.nbReadsNormal, metacluster.reads, metacluster.readsTumour, metacluster.readsNormal = metacluster.supportingReads()
            
            if 'INS' in metacluster.subclusters:
                metacluster.cv = metacluster.subclusters['INS'].cv_len()[1]

            metacluster.subclusters = None


def double_clipping_supports_INS(clusterA, clusterB, minINDELlen, technology, outDir):
    '''
    Assess if two clipping clusters support an insertion event.

    Input:
        1. clusterA: clipping cluster A
        2. clusterB: clipping cluster B
        3. minINDELlen: minimum INS and DEL lenght
        4. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        5. outDir: output file
    
    Output:
        1. boolean: clippings support an INS event (True) or not (False) 
        2. consensusFasta: fasta containing consensus read sequence spanning the insertion event. None if no insertion supporting evidences found 
    '''
    ## Search for chimeric alignments completely spanning the INS fragment
    primary, supplementary, chimeric, percSameStrand = find_chimeric_alignments(clusterA, clusterB)

    ## A) Chimeric alignment found with both pieces aligning in the same orientation
    if (primary is not None) and (percSameStrand >= 75):

        ## Search for inserted sequence at clipped clusters breakpoints
        insert = find_insertion_at_clipping_bkp(primary, supplementary)

        ## a) Inserted sequence longer than threshold
        if len(insert) >= minINDELlen:
            boolean = True            
            consensusFasta = formats.FASTA()
            consensusFasta.seqDict[primary.readName] = primary.readSeq
            
        ## b) No inserted sequence or shorter than threshold
        else:
            boolean = False
            consensusFasta = None

    ## B) Chimeric alignment NOT found -> search for complementary clippings
    else:

        ## Generate fasta files containing clusters supporting reads:
        readsA = clusterA.collect_reads() 
        readsB = clusterB.collect_reads()

        readsA_file = outDir + '/seqA.fa'
        readsB_file = outDir + '/seqB.fa'

        readsA.write(readsA_file)
        readsB.write(readsB_file)

        ## Assemble clippings based on overlap 
        contigFile = assembly.assemble_overlap(readsA_file, readsB_file, technology, outDir)
        
        ## a) Reciprocal overlap found 
        if contigFile is not None:

            boolean = True
            
            ## Read fasta with contig
            consensusFasta = formats.FASTA()
            consensusFasta.read(contigFile)             

        ## b) Reciprocal overlap not found
        else:
            boolean = False
            consensusFasta = None

    return boolean, consensusFasta


def find_chimeric_alignments(clusterA, clusterB):
    '''
    Search for chimeric read alignments connecting two clipping clusters. Select one as consensus if multiple are identified. 

    Input:
        1. clusterA: Clipping cluster object
        2. clusterB: Clipping cluster object

    Output:
        1. primary: clipping event for representative primary alignment
        2. supplementary: clipping event for representative supplementary alignment
        3. chimericSorted: list of tuples. Each tuple is composed by two clipping events corresponding to primary and supplementary alignments, respectively. 
        4. percSameStrand: percentage of chimeric alignments with both fragments having the same orientation
    ''' 
    ### 1. Identify chimeric read alignments connecting both clipping clusters
    chimeric = []
    sameStrand = 0
    oppositeStrand = 0

    # For each clipping event composing cluster A        
    for clippingA in clusterA.events:

        # For each clipping event composing cluster B
        for clippingB in clusterB.events:

            # Clipping events supported by the same read (chimeric alignments)
            if (clippingA.readName == clippingB.readName):
                
                ## Determine if fragments in chimeric alignments have the same or opposite orientations
                # A) Same
                if (clippingA.reverse == clippingB.reverse):
                    sameStrand += 1

                # B) Opposite
                else:
                    oppositeStrand += 1

                ## Determine which clipping is primary and supplementary
                # A) Clipping A primary; clipping B supplementary
                if (clippingA.supplementary == False) and (clippingB.supplementary == True):
                    primary = clippingA
                    supplementary = clippingB
                    chimeric.append((primary, supplementary))

                # B) Clipping A supplementary; clipping B primary
                elif (clippingA.supplementary == True) and (clippingB.supplementary == False):
                    primary = clippingB
                    supplementary = clippingA
                    chimeric.append((primary, supplementary))

                # C) Both clippings supplementary (Discard! These cases are not informative as supplementary alignments are hardclipped so don´t allow to pick the inserted fragment)
                
    # Exit if not chimeric alignments found
    if not chimeric:
        return None, None, None, None

    ### 2. Select the clipping events with the longest supplementary alignment as representative
    chimericSorted = sorted(chimeric, key=lambda alignments: alignments[1].refLen, reverse=True)
    primary = chimericSorted[0][0]
    supplementary = chimericSorted[0][1]

    ### 3. Compute fraction of chimeric alignments with same orientation
    percSameStrand = float(sameStrand) / (sameStrand + oppositeStrand) * 100

    return primary, supplementary, chimericSorted, percSameStrand


def find_insertion_at_clipping_bkp(primary, supplementary):
    '''
    Search for inserted DNA between the two clipping event breakpoints. 

    Input:
        1. primary: Clipping event object corresponding to the primary alignment
        2. supplementary: Clipping event object corresponding to the supplementary alignment

    Output:

        1. insert: DNA fragment inserted between clipping breakpoints
    ''' 
    ### Check if there is an unaligned piece of read sequence between the primary and supplementary alignment breakpoints
    ## A) Primary right clipped 
    # -: aligned; _: clipped            readBkpA
    # primary -----------------------------*______________________
    #                                       <<<<insert>>>>*------- supplementary
    #                                                  readBkpB == total read length -  length of the supplementary piece of sequence aligned
    if (primary.clippedSide == 'right'):
        readBkpA = primary.readBkp
        readBkpB = len(primary.readSeq) - len(supplementary.readSeq)

    ## B) Primary left clipped
    # -: aligned; _: clipped        readBkpB
    # primary        ______________________*-----------------------------
    # supplementary  -------*<<<<insert>>>>                            
    #               readBkpA == length of the supplementary piece of sequence aligned
    else:
        readBkpA = len(supplementary.readSeq)
        readBkpB = primary.readBkp             

    # Extract inserted sequence between clipping breakpoints
    insert = primary.readSeq[readBkpA:readBkpB]

    return insert


def INS_type_metaclusters(metaclusters, reference, annotations, processes, viralDb, rootOutDir):
    '''
    For each metacluster provided as input determine the type of insertion

    Input:
        1. metaclusters: list of metaclusters supporting insertion events
        2. reference: Path to the reference genome in fasta format (bwa mem and minimap2 indexes must be located in the same folder)
        3. annotations: Dictionary containing one key per type of annotation loaded and bin databases containing annotated features as values 
        4. processes: Number of processes
        5. viralDb: path to viral database in fasta format
        6. rootOutDir: Root output directory
    '''

    # Initialize variables:
    allHits_genome = {}
    groupedEntries = {}
    allHits_viral = {}

    # Dont search for viruses if there is not viralDb
    if not viralDb:
        ## 1. Create fasta containing all consensus inserted sequences 
        msg = '1. Create fasta containing all consensus inserted sequences'
        log.info(msg) 
        fastaPath = insertedSeq2fasta(metaclusters, rootOutDir)

        ## 2. Align consensus inserted sequences
        msg = '2. Align consensus inserted sequences'
        log.info(msg) 

        ## 2.1 Align consensus inserted sequences into the reference genome
        msg = '2.1 Align consensus inserted sequences into the reference genome'
        log.info(msg)    
        SAM_genome = alignment.alignment_bwa(fastaPath, reference, 'alignments_genome', processes, rootOutDir)

        ## Convert SAM to PAF
        PAF_genome = alignment.sam2paf(SAM_genome, 'alignments_genome', rootOutDir)

        ## Organize hits according to their corresponding metacluster
        allHits_genome = alignment.organize_hits_paf(PAF_genome) 

        ## 2.2 Align consensus inserted sequences into the reference genome (splicing-aware)
        msg = '2.2 Align consensus inserted sequences into the reference genome (splicing-aware)'
        log.info(msg)    

        ## Minimap index for the reference
        index = os.path.splitext(reference)[0] + '.mmi'
        SAM_splicing = alignment.alignment_minimap2_spliced(fastaPath, index, 'alignments_spliced', processes, rootOutDir)

        ## Convert SAM to BAM
        BAM_splicing = bamtools.SAM2BAM(SAM_splicing, rootOutDir)

        ## Convert BAM to BED
        BED_path = bamtools.BAM2BED(BAM_splicing, rootOutDir)

        ## Organize hits according to their corresponding metacluster
        allHits_splicing = formats.BED()
        allHits_splicing.read(BED_path, 'List', None)
        groupedEntries = allHits_splicing.group_entries_by_name()

    # Search only for viruses if there is viralDb
    else:
        ## 1. Create fasta containing all consensus inserted sequences and checking sequences local complexity
        msg = '1. Create fasta containing all consensus inserted sequences'
        log.info(msg) 
        fastaPath = insertedSeq2fasta(metaclusters, rootOutDir, lccFilter=True)

        ## 2. Align consensus inserted sequences
        msg = '2. Align consensus inserted sequences'
        log.info(msg) 

        ## 2.3 Align consensus inserted sequences into the viral database
        msg = '2.1 Align consensus inserted sequences into the viral database'
        log.info(msg)  

        #start_time = time.time()
        SAM_viral = alignment.alignment_bwa(fastaPath, viralDb, 'alignments_viral', processes, rootOutDir)
        #print("--- %s seconds SAM_viral ---" % (time.time() - start_time))
        
        ## Convert SAM to PAF
        PAF_viral = alignment.sam2paf(SAM_viral, 'alignments_viral', rootOutDir)

        ## Organize hits according to their corresponding metacluster
        allHits_viral = alignment.organize_hits_paf(PAF_viral)
    
    ## 3. For each metacluster determine the insertion type
    msg = '3. For each metacluster determine the insertion type'
    log.info(msg)   
    
    # For each metacluster
    for metacluster in metaclusters:

        ## 3.1 Collect consensus inserted sequence hits
        metaId = str(metacluster.ref) + ':' + str(metacluster.beg) + '-' + str(metacluster.end)

        ## Hits in the reference genome
        if metaId in allHits_genome:
            hits_genome = allHits_genome[metaId]
        
        else:
            hits_genome = formats.PAF()

        ## Hits in the reference genome (splice-aware alignment)
        if metaId in groupedEntries:
            hits_splicing = groupedEntries[metaId]

        else:
            hits_splicing = []

        ## Hits in the viral database
        if metaId in allHits_viral:
            hits_viral = allHits_viral[metaId]

        else:
            hits_viral = formats.PAF()

        ## 3.2 Insertion type inference
        if not viralDb:
            metacluster.determine_INS_type(hits_genome, hits_splicing, hits_viral, annotations['REPEATS'], annotations['TRANSDUCTIONS'], annotations['EXONS'])
        else:
            # Dont look if it was already identified
            # TODO: This has to be organised in another way.
            if 'IDENTITY' not in metacluster.SV_features.keys():
                metacluster.determine_INS_type(hits_genome, hits_splicing, hits_viral, annotations['REPEATS'], annotations['TRANSDUCTIONS'], annotations['EXONS'], types2Search=['VIRUS'])


def structure_inference_parallel(metaclusters, consensusPath, transducedPath, transductionSearch, processes, rootDir):
    '''
    Infer structure for a list of INS metacluster objects. Parallelize by distributing metaclusters by processes. 

    Input:
        1. metaclusters: list of metacluster objects
        2. consensusPath: path to fasta file containing retrotransposon consensus sequences
        3. transducedPath: path to fasta containing transduced sequences downstream of source elements
        4. transductionSearch: boolean specifying if transduction search is enabled (True) or not (False)
        5. processes: number of processes
        6. rootDir: Root output directory
    
    Output:
        1. metaclustersOut: list of metacluster objects with structure information stored at 'SV_features' dict attribute
    '''
    ## 1. Create tuple list for multiprocessing
    msg = '1. Create tuple list for multiprocessing'
    log.subHeader(msg)      
    metaclustersOut = []
    tupleList = []

    for metacluster in metaclusters:
        
        ## Skip structure inference if insertion type not available or not solo, partnered or orphan transduction
        # Note: investigate why INS_TYPE is not defined in some metaclusters
        if ('INS_TYPE' not in metacluster.SV_features) or (metacluster.SV_features['INS_TYPE'] not in ['solo', 'partnered', 'orphan']):
            metaclustersOut.append(metacluster) 
            continue

        ## Add to the list
        fields = (metacluster, consensusPath, transducedPath, transductionSearch, rootDir)
        tupleList.append(fields)

    ## 2. Infer structure
    msg = '2. Infer structure'
    log.subHeader(msg)       

    pool = mp.Pool(processes=processes)
    metaclustersStructure = pool.starmap(structure_inference, tupleList)
    pool.close()
    pool.join()
      
    ## 3. Create final list with metaclusters to be reported as output
    msg = '3. Create final list with metaclusters to be reported as output'
    metaclustersOut = metaclustersOut + metaclustersStructure 

    return metaclustersOut

def structure_inference(metacluster, consensusPath, transducedPath, transductionSearch, rootDir):
    '''
    Wrapper to call 'determine_INS_structure' method for a given INS metacluster provided as input

    Input:
        1. metacluster: INS metacluster 
        2. consensusPath: path to fasta file containing retrotransposon consensus sequences
        3. transducedPath: path to fasta containing transduced sequences downstream of source elements
        4. transductionSearch: boolean specifying if transduction search is enabled (True) or not (False)
        5. rootDir: Root output directory
    
    Output:
        1. metacluster: INS metacluster containing structural properties  
    '''
    # Create output directory
    metaInterval = '_'.join([str(metacluster.ref), str(metacluster.beg), str(metacluster.end)])
    outDir = rootDir + '/' + metaInterval
    unix.mkdir(outDir)

    # Infer structure
    structure = metacluster.determine_INS_structure(consensusPath, transducedPath, transductionSearch, outDir)

    # Add structure info to the metacluster
    metacluster.SV_features.update(structure) 
    
    # Remove output directory
    unix.rm([outDir])

    return metacluster

def insertedSeq2fasta(metaclusters, outDir, lccFilter=False):
    '''
    Collect all the consensus inserted sequences from a list of metaclusters supporting INS and 
    generate fasta file containing them 

    Input:
        1. metaclusters: list of metaclusters
        2. outDir: output directory
        3. lccFilter: Boolean. Filter reads by local complexity. Default= False.

    Output:
        1. fastaPath: fasta file containing all the consensus inserted sequences
    '''
    ## Initiate FASTA object
    FASTA = formats.FASTA()

    # For each metacluster
    for metacluster in metaclusters:

        ## Skip insertion type inference if consensus event not available
        if metacluster.consensusEvent is None:
            continue               

        ## Retrieve inserted sequence and add to the FASTA
        metaclusterId = metacluster.ref + ':' + str(metacluster.beg) + '-' + str(metacluster.end)
        insert = metacluster.consensusEvent.pick_insert()

        # NOTE 2020: New 2020
        if lccFilter:
            complexity = Bio.SeqUtils.lcc.lcc_simp(insert)
            print ('complexity ' + str(insert) + ' '+ str(complexity))

            # NOTE 2020: If works, put as option
            if complexity > 1.9:
                FASTA.seqDict[metaclusterId] = insert
        
        FASTA.seqDict[metaclusterId] = insert
        
    ## Write fasta         
    fastaPath = outDir + '/inserted_sequences.fa'
    FASTA.write(fastaPath)    

    return fastaPath


def assignAligments2metaclusters_paf(metaclusters, PAF_path):
    '''
    Map alignments to their corresponding metacluster. 

    Input:
        1. metaclusters: list of metaclusters
        2. PAF_path: Path to path file containing alignments to asign

    Output:
        1. tupleList: List of tuples. Each tuple contain two elements: 
            1) metacluster object  
            2) PAF object containing the corresponding alignments for the metacluster consensus inserted sequence
    '''
    ## 1. Create a dictionary to organize the data
    hits = {}

    for metacluster in metaclusters:
        metaclusterId = metacluster.ref + ':' + str(metacluster.beg) + '-' + str(metacluster.end)
        hits[metaclusterId] = (metacluster, formats.PAF())

    ## 2. Read PAF file and add hits to the metaclusters
    ## Read PAF 
    PAF = formats.PAF()
    PAF.read(PAF_path)

    # For each read alignment 
    for alignment in PAF.alignments:
        hits[alignment.qName][1].alignments.append(alignment)
    
    ## 3. Generate list of tuples
    tupleList = list(hits.values())

    return tupleList


def assignAligments2metaclusters_sam(metaclusters, SAM_path):
    '''
    Map alignments to their corresponding metacluster. 

    Input:
        1. metaclusters: list of metaclusters
        2. SAM_path: Path to SAM file containing alignments to asign

    Output:
        1. metaclustersHits: Nested list. Each list element is composed by a list with two elements: 
            1) metacluster object  
            2) List of aligned segment objects
    '''
    ## 1. Create a dictionary to organize the data
    hits = {}

    for metacluster in metaclusters:
        metaclusterId = metacluster.ref + ':' + str(metacluster.beg) + '-' + str(metacluster.end)
        hits[metaclusterId] = [metacluster, []]

    ## 2. Read SAM file and add hits to the metaclusters
    ## Read SAM 
    SAM = pysam.AlignmentFile(SAM_path, "r")

    # For each read alignment 
    for alignment in SAM:
        hits[alignment.query_name][1].append(alignment)
    
    ## 3. Generate nested list
    metaclustersHits = list(hits.values())

    return metaclustersHits

def search4bridges_metaclusters_parallel(metaclusters, maxBridgeLen, minMatchPerc, minReads, minPercReads, annotations, refDir, viralDb, processes, rootDir):
    '''
    Search for transduction or repeat bridges at BND junctions for a list of metacluster objects

    Input:
        1. metaclusters: list of input metacluster objects supporting BND
        2. maxBridgeLen: maximum supplementary cluster length to search for a bridge
        3. minMatchPerc: minimum percentage of the supplementary cluster interval to match in a transduction or repeats database to make a bridge call
        4. minReads: minimum number of reads supporting the bridge
        5. minPercReads: minimum percentage of clipping cluster supporting reads composing the bridge
        6. annotations: Dictionary containing one key per type of annotation loaded and bin databases containing annotated features as values 
        7. refDir: directory containing reference databases. 
        8. processes: number of processes
        9. rootDir: root output directory

    For each metacluster update 'bridgeClusters' and 'bridgeType' attributes
    '''    

    if not viralDb:
        ## 1. Generate index containing consensus retrotranposon sequences + source elements downstream regions
        ## Consensus retrotransposon sequences
        consensusPath = refDir + '/consensusDb.fa'
        consensus = formats.FASTA()
        consensus.read(consensusPath)

        ## Transduced regions
        transducedPath = refDir + '/transducedDb.fa.masked'
        transduced = formats.FASTA()
        transduced.read(transducedPath)
        
        ## Merge in a single fasta
        allSeqs = formats.FASTA()
        allSeqs.seqDict = {**consensus.seqDict, **transduced.seqDict}

        ## Write fasta
        fastaPath = rootDir + '/reference_sequences.fa'
        allSeqs.write(fastaPath)

        ## Index fasta 
        fileName = 'reference_sequences'  
        index = alignment.index_minimap2(fastaPath, fileName, rootDir)

    # When analysing viruses no index is needed:
    else:
        index = None

    ## 2. Create tuple list for multiprocessing
    tupleList = []

    for metacluster in metaclusters:
        
        ## Add to the list
        fields = (metacluster, maxBridgeLen, minMatchPerc, minReads, minPercReads, annotations, index, viralDb, rootDir)
        tupleList.append(fields)

    ## 3. Search for bridges
    pool = mp.Pool(processes=processes)
    metaclusters = pool.starmap(search4bridges_metacluster, tupleList)
    pool.close()
    pool.join()

    return metaclusters

def search4bridges_metacluster(metacluster, maxBridgeLen, minMatchPerc, minReads, minPercReads, annotations, index, viralDb, rootDir):
    '''
    Search for transduction or repeat bridges at BND junctions for a metacluster object

    Input:
        1. metacluster: metacluster object
        2. maxBridgeLen: maximum supplementary cluster length to search for a bridge
        3. minMatchPerc: minimum percentage of the supplementary cluster interval to match in a transduction or repeats database to make a bridge call
        4. minReads: minimum number of reads supporting the bridge
        5. minPercReads: minimum percentage of clipping cluster supporting reads composing the bridge
        6. annotations: dictionary containing one key per type of annotation loaded and bin databases containing annotated features as values (None for those annotations not loaded)
        7. index: minimap2 index for consensus retrotransposon sequences + transduced regions database
        8. rootDir: root output directory

    Output:
        1. metacluster: metacluster object with updated bridge information
    '''
    ## 1. Create output directory
    metaInterval = '_'.join([str(metacluster.ref), str(metacluster.beg), str(metacluster.end)])
    outDir = rootDir + '/' + metaInterval
    unix.mkdir(outDir)

    ## 2. Search for bridge
    metacluster.search4bridge(maxBridgeLen, minMatchPerc, minReads, minPercReads, annotations, index, viralDb, outDir)

    ## 3. Remove output directory        
    unix.rm([outDir])

    return metacluster


def search4junctions_metaclusters(metaclusters, refLengths, processes, minReads, minPercReads, reference, refDir, viralDb, outDir):
    '''
    Search for BND junctions for a list of input metacluster objects of the type BND.

    Each BND junction will correspond to a connection between two BND metaclusters
    
    BND junctions are identified based on the analysis of supplementary alignments.

    Input:
        1. metaclusters: list of input metacluster objects supporting BND
        2. refLengths: dictionary containing reference ids as keys and as values the length for each reference. 
        3. processes: number of processes
        4. minReads: minimum number of reads supporting the BND junction 
        5. minPercReads: minimum percentage of metacluster and partner metacluster reads supporting the BND junction
        6. refDir
        7. outDir: output directory

    Output:
        1. allJunctions: list of BND junction objects
    '''
    ## 1. Organize metaclusters into a dictionary
    metaclustersDict = events.events2nestedDict(metaclusters, 'BND')

    ## 2. Organize metaclusters into a bin database
    metaclustersBinDb = structures.create_bin_database_parallel(refLengths, metaclustersDict, processes)

    ## 3. Search for BND junctions for each metacluster
    allJunctions = []
    includedJunctions = []

    ## For each metacluster
    for metacluster in metaclusters:
    
        # Search for junctions
        junctions = metacluster.search4junctions(metaclustersBinDb, minReads, minPercReads)
        
        # Add junctions to the final list avoiding redundancies
        for junction in junctions:

            # Skip redundant junctions 
            if (junction.junctionCoord() in includedJunctions) or (junction.junctionCoord_rev() in includedJunctions):
                continue

            # Add junction to the list
            allJunctions.append(junction)
            includedJunctions.append(junction.junctionCoord())
    
    return allJunctions

def analyse_BNDjunction_structure(junction, reference, viralDb, outDir):
    '''
    Function to analyse BND_junction structure
    Input:
        - allJunctions: List containing BND_junctions objects
    Output:
        - This function has no output. It changes junction atributes:
            - junction.junctionConsSeq
            - junction.bridge.bridgeType
            - junction.bridge.family
            - junction.bridge.srcId
            - junction.bkpsConsSupport
    '''

    junctionOutDir = outDir +'/'+ junction.junctionCoord() 
    unix.mkdir(junctionOutDir)

    junctionPAFChain = None

    #for junction in allJunctions:
    # TODO: REMOVE this line: junction = BND_junction
    # TODO: REMOVE this line: junction.bridge = <clusters.BRIDGE object at 0x7f7d9cf8f940> antes: {'TRANSDUCTION': {'6p24.1': <clusters.BRIDGE object at 0x7f891bff6b38>}, 'REPEAT': {}}

    # Pick reads supporting the BND_junction  
    junctionsList = junction.extractSupportingRead()

    # Pick the event that is present in both metaclusters and all bridges (which is the one with the highest number in junctionsList[1]) and, from those, pick the one with the lowest number os suppAlignments (which is the one with the lowest number in junctionsList[2])
    # If there are many that match these conditions, just pick the first one.
    # Pick metacluster events of the read that will be used as a template
    metaclustersEvents = junctionsList[0][0]

    # Pick the read coordinates that delimit the BND_junction
    # TODO: int instead of a list
    junctionInterval=[]
    # Pick event which sequence will be used as template
    templateEvent = junctionsList[0][0][0]
    junctionInterval.append(templateEvent.readBkp)

    # Get consensus sequence from all reads
    polishedFastaEntireSequence = get_consensus_BNDjunction(templateEvent, junctionsList, junctionOutDir)

    # TODO: remove junctionsList
    # Get consensus sequence that spans the BND_junction +- 1000bp
    # NOTE 2020: New EGA 2020
    # TODO 2020: Put these checks in a good way
    if polishedFastaEntireSequence:
        if os.stat(polishedFastaEntireSequence).st_size !=0 :
            polishedFastaInterval, polishedFastaIntervalObj = junction.get_consensus_BNDInterval(polishedFastaEntireSequence, junctionInterval, junctionOutDir)

            # Get reference coordinates to perform target alignment
            targetIntervalList = []
            # metaclusterA
            offset = 1000
            intervalBeg = junction.metaclusterA.bkpPos - offset if junction.metaclusterA.bkpPos - offset >= 0 else 0 ## Set lower bound
            intervalEnd = junction.metaclusterA.bkpPos + offset if junction.metaclusterA.bkpPos + offset >= 0 else 0 ## Set lower bound
            intervalCoord = junction.metaclusterA.ref + ':' + str(min([intervalBeg,intervalEnd])) + '-' + str(max([intervalBeg,intervalEnd]))

            targetIntervalList.append(intervalCoord)
            # metaclusterB
            # TODO: PUT THIS IN A FUNCTION SO WE DONT HAVE TO REPEAT IT
            offset = 1000
            intervalBeg = junction.metaclusterB.bkpPos - offset if junction.metaclusterB.bkpPos - offset >= 0 else 0 ## Set lower bound
            intervalEnd = junction.metaclusterB.bkpPos + offset if junction.metaclusterB.bkpPos + offset >= 0 else 0 ## Set lower bound
            intervalCoord = junction.metaclusterB.ref + ':' + str(min([intervalBeg,intervalEnd])) + '-' + str(max([intervalBeg,intervalEnd]))
            

            targetIntervalList.append(intervalCoord)

            # TODO: fix reference file
            #target = sequences.create_targeted_fasta(targetIntervalList, reference, outDir)
            # Target reference from coordinates to perform target alignment

            target = sequences.create_targeted_fasta(targetIntervalList, reference, junctionOutDir)

            junctionPAFChain = junction.get_BNDjunction_chain(polishedFastaInterval, target, viralDb, junctionOutDir)

    if junctionPAFChain:
        aligTNames, aligCoordinates = junction.get_consensus_BNDSeq(junctionPAFChain, polishedFastaIntervalObj)
        junction.analyse_BNDjunction_bridge(aligTNames, aligCoordinates, polishedFastaIntervalObj, viralDb, junctionOutDir)
        junction.consensusBNDjunction_support_Bkps(junctionPAFChain)

        '''
        # TODO: think if necessary
        # Check if the pieces are close or there is space between them
        alignmentsGaps = []
        for alig in junctionPAFChain.alignments:
            beg = alig.qBeg
            if beg == junctionPAFChain.interval()[0]:
                end = alig.qEnd
            else:
                beg = alig.qBeg
                alignmentsGaps.append(beg-end)
                end = alig.qEnd

        if max(alignmentsGaps) > 100:
            print ('Gap')
        else:
            print ('No gap')
        '''
        # TODO: Check if bridge element correponds to bridge type: if not WARNING.
    #unix.rm([junctionOutDir])

    return junction

def get_consensus_BNDjunction(templateEvent, junctionsList, outDir):
    '''
    Function to get consensus sequence that completely spans a BND_junction
    Input:
        - templateEvent: event supporting the BND_junction that will be used as template
        - junctionsList:
                - field 0: list of metaclusters that contain the same read
                - field 1: times that the read appears in a different events
                - field 2: number of supplementary alignments of the read.
        - outDir
    Output:
        - Fasta file containing consensus sequence that completely spans a BND_junction
    '''

    # 1. Make fasta with template sequence
    templateFastaObj = formats.FASTA()
    templateFastaObj.seqDict[templateEvent.readName] = templateEvent.readSeq
    templateFastaPath = outDir + '/template_sequence.fa'
    templateFastaObj.write(templateFastaPath)

    # 2. Make fasta with supporting reads
    supportingReadsFastaObj = formats.FASTA()
    for sublist in junctionsList[1:]:
        supportingReadsFastaObj.seqDict[sublist[0][0].readName] = sublist[0][0].readSeq
    
    supportingReadsFastaPath = outDir + '/supporting_sequences.fa'
    supportingReadsFastaObj.write(supportingReadsFastaPath)

    # 3. POlish with racon
    # TODO: change the following lines: 
    #polishedFasta = assembly.polish_racon(templateFastaPath, supportingReadsFastaPath, confDict['technology'], confDict['rounds'], outDir)
    polishedFastaEntireSequence = assembly.polish_racon(templateFastaPath, supportingReadsFastaPath, 'NANOPORE', 1, outDir)

    return polishedFastaEntireSequence
    
## CLASSES ##
    

#############
## CLASSES ##
#############

class cluster():
    '''
    Events cluster class. A cluster is composed by a set of events.
    Each event is supported by a single read. One cluster can completely represent a single structural
    variation event or partially if multiple clusters are required (see 'metaCluster' class)
    '''
    number = 0 # Number of instances

    def __init__(self, events, clusterType):
        '''
        '''
        cluster.number += 1 # Update instances counter
        self.id = 'CLUSTER_' + str(cluster.number)
        self.clusterId = None

        # Define list of events composing the cluster and cluster type
        self.events = events
        self.clusterType = clusterType

        # Set cluster's reference, begin and end position
        self.ref, self.beg, self.end = self.coordinates() 

        # Cluster filtering
        self.filters = None
        self.failedFilters = None
        self.nbOutliers = 0

        # Update event's clusterId attribute
        for event in self.events:
            event.clusterId = self.id        

    def length(self):
        '''
        Compute cluster interval length
        '''
        return self.end - self.beg
        
    def sort(self):
        '''
        Sort events in increasing coordinates order
        '''
        self.events.sort(key=lambda event: event.beg)

    def sort_by_length(self):
        '''
        Sort events in increasing length ordering
        '''
        return sorted(self.events, key=lambda event: event.length)
        
    def coordinates(self):
        '''
        Compute cluster ref, beg and end coordinates. 
        
        Begin and end will correspond to the left and rightmost positions, respectively
        '''
        # Sort events from lower to upper beg coordinates
        self.sort()  
        
        # Define cluster coordinates 
        ref = self.events[0].ref
        beg = min([event.beg for event in self.events])
        end = max([event.end for event in self.events])
            
        return ref, beg, end

    def add(self, events2add):
        '''
        Incorporate events into the cluster and redefine cluster beg and end
        positions accordingly

        Input:
            1. events2add: List of events to be added to the cluster
        '''
        # Add events to the cluster  
        self.events = self.events + events2add

        # Resort and redefine cluster begin and end coordinates
        self.ref, self.beg, self.end = self.coordinates() 

        # Update event's clusterId attribute
        for event in events2add:
            event.clusterId = self.id

    def remove(self, events2remove):
        '''
        Remove list of events from the cluster and redefine cluster beg and end
        positions accordingly

        Input:
            1. events2remove: List of events to be removed 
        '''
        ## 1. Remove events from the metacluster ##
        self.events = [event for event in self.events if event not in events2remove]

        ## 2. Resort and redefine metacluster begin and end coordinates ##
        if self.events:
            self.ref, self.beg, self.end = self.coordinates()

    def pick_median_length(self):
        '''
        Return event whose length is at the median amongst all cluster supporting events
        '''
        ## Sort events by their length
        sortedEvents = self.sort_by_length()

        ## Compute the index for the event with the median length
        median = (len(sortedEvents) - 1)/2  # minus 1 because the first element is index 0
        medianIndex = int(math.ceil(median))

        ## Pick event located at the median 
        return sortedEvents[medianIndex]

    def cv_len(self):
        '''
        Compute mean length and coefficient of variation (cv)

        Output:
            1. meanLen: Mean length
            2. cv: Coefficient of variation 
        '''
        lengths = [ event.length for event in self.events]
        meanLen = np.mean(lengths)
        std = np.std(lengths)
        cv = std / meanLen * 100 

        return meanLen, cv

    def collect_reads(self):
        '''
        Create FASTA object containing cluster supporting reads.
        
        Output:
            1. FASTA: FASTA object containing cluster supporting reads
        '''
        ## Initiate FASTA object
        FASTA = formats.FASTA()

        ## Add reads supporting the events to the FASTA
        for event in self.events:

            FASTA.seqDict[event.readName] = event.readSeq

        return FASTA
    
    def nbReads(self):
        '''
        Return the total number of reads supporting the cluster and the list of read ids
        '''
        readList = list(set([event.readName for event in self.events]))
        nbReads = len(readList)
        return nbReads, readList

    def supportingReads(self):
        '''
        Compute the total number of cluster supporting reads and generate a list of supporting read ids

        Output:
            1. nbTotal: Total number of cluster supporting reads
            2. nbTumour: Number of cluster supporting reads in the tumour 
            3. nbNormal: Number of cluster supporting reads in the normal
            4. reads: List containing all the cluster supporting reads
            5. readsTumour: List containing cluster supporting reads in the tumour sample
            6. readsNormal: List containing cluster supporting reads in the normal sample
        '''
        reads = []
        readsTumour = []
        readsNormal = []

        ## 1. Create non-redundant list of metacluster supporting reads
        for event in self.events: 

            # Add read name not already included in the list  
            if event.readName not in reads:
                reads.append(event.readName)

            ## Tumour and matched normal
            # a) Event identified in the TUMOUR sample
            if event.sample == "TUMOUR":
            
                # Add read name not already included in the list  
                if event.readName not in readsTumour:
                    readsTumour.append(event.readName)
            
            # b) Event identified in the matched NORMAL sample
            elif event.sample == "NORMAL":

                # Add read name not already included in the list  
                if event.readName not in readsNormal:
                    readsNormal.append(event.readName)            

        ## 2. Compute number of supporting reads
        nbTotal = len(reads) # total

        # a) Unpaired mode
        if self.events[0].sample is None:
            nbTumour, nbNormal, readsTumour, readsNormal = [None, None, None, None]

        # b) Paired mode
        else:

            nbTumour = len(readsTumour) # tumour
            nbNormal = len(readsNormal) # normal

        return nbTotal, nbTumour, nbNormal, reads, readsTumour, readsNormal
        
    def nbEvents(self):
        '''
        Return the number of events composing the cluster
        '''
        nbTumour = 0
        nbNormal = 0

        for event in self.events:
            # a) Event identified in the TUMOUR sample
            if event.sample == "TUMOUR":
                nbTumour += 1
            
            # b) Event identified in the matched NORMAL sample
            elif event.sample == "NORMAL":
                nbNormal += 1
            
            # c) SINGLE sample mode
            else:
                nbTumour = None
                nbNormal = None
                break

        nbTotal = len(self.events)

        return nbTotal, nbTumour, nbNormal
    
    def dupPercentage(self):
        '''
        Return the percentage of duplicates in cluster
        '''
        
        nbDup = 0
        nbTotal = len(self.events)
        
        for event in self.events:
            
            if event.isDup == True:
                nbDup += 1
        
        dupPercentage = stats.fraction(nbDup, nbTotal) * 100
        
        return dupPercentage
        


class INS_cluster(cluster):
    '''
    Insertion (INS) cluster subclass
    '''
    def __init__(self, events):

        cluster.__init__(self, events, 'INS')

    def correct_fragmentation(self):
        '''
        Correct fragmented alignments over INS
        
        Before merging:
        ############<<<INS>>>##<<INS>>###<<<INS>>>##############
        
        After merging:
        ############<<<<<<<<INS>>>>>>>>##############        
        '''
        ## 1. Organize INS events into a dictionary according to their supporting read
        eventsByReads =  {}

        for INS in self.events:
        
            # a) First event supported by that read
            if INS.readName not in eventsByReads:
                eventsByReads[INS.readName] = [INS]

            # b) There are previous events supported by that read
            else:
                eventsByReads[INS.readName].append(INS)
        
        ## 2. Merge INS events supported by the same read
        merged_list = []
        fragmented_list = []

        # For each read
        for readId, INS_list in eventsByReads.items():
            
            ## Read supporting multiple INS events
            if len(INS_list) > 1:

                ## 2.1 Do merging of fragmented INS
                merged = events.merge_INS(INS_list)

                ## 2.2 Add merged INS
                merged_list.append(merged)

                ## 2.3 Update fragmented alignments list
                fragmented_list = fragmented_list + INS_list
    
        ## 3. Update cluster
        self.add(merged_list)
        self.remove(fragmented_list)

    def identify_subclusters(self, minClusterSize):
        '''
        Identify subclusters of INS with similar lengths

        Input:
            1. minClusterSize: minimum size to create a subcluster

        Output:
            1. subclusters: list of subclusters
        '''

        subclusters = []

        ## Compute metrics based on events length
        meanLen, cv = self.cv_len()

        ## A) Cluster with heterogeneous lengths  
        if cv > 15:

            ## 1. Cluster events according to their length
            clustered = clustering.KMeans_clustering(self.events, None, 'length')

            ## 2. Filter out subclusters supported by less than X events
            clusteredFiltered = []

            for cluster in clustered.values():

                if len(cluster) >= minClusterSize:
                    clusteredFiltered.append(cluster)

            ## 3. Create subclusters
            # NOTE: Don´t create subclusters if:
            #   a) Cluster no fragmented into any subcluster passing support filter OR
            #   b) Cluster fragmented in more than 2 subclusters passing support filter (indicative of noisy region)
            if len(clusteredFiltered) in [1, 2]:

                for events in clusteredFiltered:

                    subcluster = INS_cluster(events)
                    subclusters.append(subcluster)

        return subclusters     


class DEL_cluster(cluster):
    '''
    Deletion (DEL) cluster subclass
    '''
    def __init__(self, events):

        cluster.__init__(self, events, 'DEL')


class CLIPPING_cluster(cluster):
    '''
    Clipping cluster subclass
    '''
    def __init__(self, events):

        cluster.__init__(self, events, 'CLIPPING')
        self.supplCluster = None
        self.orientation = None  

    def infer_breakpoint(self):
        '''
        Infer clipping cluster breakpoint position

        Output:
            1. bkpPos: consensus breakpoint position
        '''
        ## Retrieve all positions
        positions = [clipping.beg for clipping in self.events]

        ## Count the number of reads supporting each position
        positionCounts = Counter(positions)

        ## Organize positions by degree of support into a dictionary
        supportDict = {}

        for position in positionCounts:
            nbReads = positionCounts[position]

            if nbReads not in supportDict:
                supportDict[nbReads] = [position]
            else:
                supportDict[nbReads].append(position)
        
        ## Sort by level of support
        supportLevels = sorted(supportDict.keys(), reverse=True)

        ## Select breakpoints with maximum support level
        maxSupportLevel = supportLevels[0]
        maxSupportBkps = supportDict[maxSupportLevel]
        nbBkps = len(maxSupportBkps)

        ## Resolve breakpoint
        # A) No ambiguity (single breakpoint maximum read support)
        if nbBkps == 1:
            bkpPos = maxSupportBkps[0]

        # B) Ambiguous bkp (several breakpoints with maximum read support)
        else:
            ## Compute median
            medianPos = int(np.median(positions))

            ## Compute breakpoint distance to the median
            dist2medians = []

            for bkp in maxSupportBkps:

                dist2median = abs(bkp - medianPos)
                dist2medians.append((dist2median, bkp))

            ## Sort breakpoints in decreasing distance order
            dist2medians = sorted(dist2medians, key=itemgetter(0))

            ## Select breakpoint closest to the median 
            # Note: use one arbitrary if still several possibilities
            bkpPos = dist2medians[0][1]

        return bkpPos

    def search4complAlignments(self):
        '''
        Search for complementary alignments for each clipping event composing the clipping cluster.
        
        A complementary alignment span a piece of the clipped read sequence, so contributes to explain 
        the full conformation of the read alignment. 

        Output:
            1. allComplementaryAlignments: complementary alignments identified for all the clipping events composing the cluster. 
            These are organized into a dictionary with the references as keys and the list of complementary alignments as values 
        '''
        ## 1. Collect complementary alignments for each clipping event            
        complementaryAlignmentsDictList = []

        for clipping in self.events:

            ## 1.1 Retrieve suppl. alignments            
            supplAlignmentsDict = clipping.parse_supplAlignments_field()

            ## 1.2 Organize suppl. alignments into a simple list
            supplAlignments = structures.dict2list(supplAlignmentsDict)

            ## 1.3 Search for complementary alignments
            complementaryAlignments = clipping.search4complAlignments(supplAlignments)

            ## 1.4 Organize complementary alignments into a dictionary
            complementaryAlignmentsDict = events.events2Dict(complementaryAlignments) 
            complementaryAlignmentsDictList.append(complementaryAlignmentsDict)

        ## 2. Merge complementary alignments into a single dictionary
        allComplementaryAlignments = structures.merge_dictionaries(complementaryAlignmentsDictList)

        return allComplementaryAlignments

    def cluster_suppl_positions(self, supplAlignmentsDict):
        '''
        Cluster supplementary alignments based on their begin/end alignment positions

        Input:
            1. supplAlignmentsDict: dictionary with the references as keys and the list of supplementary alignments as values 

        Output:
            1. supplClusters: Dictionary containing references as keys and nested lists of supplementary 
            alignment clusters as values.
        '''
        supplClusters = {}

        # For each reference
        for ref in supplAlignmentsDict:

            # Collect suppl. alignments list
            supplAlignments = supplAlignmentsDict[ref]

            ## Cluster suppl. alignments based on their beg and end alignment positions
            clustersBeg = clustering.distance_clustering_targetPos(supplAlignments, 100, 'beg')
            clustersEnd = clustering.distance_clustering_targetPos(supplAlignments, 100, 'end')

            ## Determine bkp side based on the biggest cluster
            biggestLenBeg = max([len(cluster.events) for cluster in clustersBeg])
            biggestLenEnd = max([len(cluster.events) for cluster in clustersEnd])

            # a) Bkp at the beg of supplementary alignment interval
            if biggestLenBeg >= biggestLenEnd:
                clusters = clustersBeg
            
            # b) Bkp at the end of supplementary alignment interval
            else:
                clusters = clustersEnd
                
            ## Add clusters to the dictionary
            supplClusters[ref] = clusters

        return supplClusters

class SUPPLEMENTARY_cluster(cluster):
    '''
    Supplementary alignment cluster subclass
    '''
    def __init__(self, events):

        cluster.__init__(self, events, 'SUPPLEMENTARY')
        self.representative = None
        self.bkpSide = None
        self.annot = None
        self.bridge = False
        self.bridgeInfo = {}
        self.orientation = None
    
    def clippOrientation(self):
        '''
        Determine cluster bkp side
        '''
        # 1. Collect supplementary events clipping sides
        [event.clippingSide() for event in self.events]
        clipSides = set([event.clipSide for event in self.events])
        
        # 2. Set supplementary cluster bkp side
        if clipSides == {'left', 'right'}:
            self.orientation = 'RECIPROCAL'
            
        elif clipSides == {'right'}:
            self.orientation = 'PLUS'
            
        elif clipSides == {'left'}:
            self.orientation = 'MINUS'

    def inferBkp_shortReads(self):
        '''
        Compute and return breakpoint position for PE short reads
        '''
        # Determine cluster bkp side
        if self.orientation == None:
            self.clippOrientation()
        
        # a) Breakpoint on the left
        if self.orientation == 'MINUS':
            coordList = [event.beg for event in self.events]
            bkpPos = max(set(coordList), key=coordList.count)

        # b) Breakpoint on the right
        elif self.orientation == 'PLUS':
            coordList = [event.end for event in self.events]
            bkpPos = max(set(coordList), key=coordList.count)
        
        # c) If orientation reciprocal or unknown, return coordinate with max nb of counts
        else:
            coordList = [event.beg for event in self.events] + [event.end for event in self.events]
            bkpPos = max(set(coordList), key=coordList.count)
        
        return bkpPos

    def bkpPos(self):
        '''
        Compute and return breakpoint position
        '''
        # a) Breakpoint on the left
        if self.bkpSide == 'beg':
            bkpPos = self.beg

        # b) Breakpoint on the right
        else:
            bkpPos = self.end
        
        return bkpPos
        
    def clusterIndex(self):
        '''
        Compute supplementary cluster read level index. Index relative to the clipping
        cluster

        ########################################## Hipothetical consensus read
        ------------- ------------ ---------------
         -----------   ----------   -------------
        ------------   ----------   -------------
      clipping_cluster suppl_cluster suppl_cluster
                        (index: 0)    (index: 1)

        Output:
            1. index: read level index relative to the clipping cluster
        '''
        ## 1. Generate list containing index position for all the suppl. alignments composing the cluster
        indexes = [supplementary.readIndex for supplementary in self.events]

        ## 2. Count number occurrences for each index
        occurences = [[index ,indexes.count(index)] for index in set(indexes)]

        ## 3. Sort indexes in decreasing order of occurrences
        occurences.sort(key = lambda x: x[1], reverse = True) 

        # 4. Select index with maximum number of instances
        index = occurences[0][0]

        ## Note: return an ambiguous flag if several possible maximum
        return index

    def support_bridge(self, maxBridgeLen, minMatchPerc, minReads, minPercReads, annotations, index, viralDb, outDir):
        '''
        Assess if supplementary cluster supports bridge or not. 

        Input:
            1. maxBridgeLen: maximum supplementary cluster length to search for a bridge
            2. minMatchPerc: minimum percentage of the supplementary cluster interval to match in a transduction or repeats database to make a bridge call
            3. minReads: minimum number of reads supporting the bridge
            4. minPercReads: minimum percentage of clipping cluster supporting reads composing the bridge
            5. annotations: dictionary containing one key per type of annotation loaded and bin databases containing annotated features as values (None for those annotations not loaded)
            6. index: minimap2 index for consensus retrotransposon sequences + source element downstream regions
            7. outDir: output directory

        Output: Update 'bridge' and 'bridgeInfo' attributes
        '''

        ## 1. Search for unaligned bridge sequence at BND junction (algorithm 1)
        self.bridge, supportType, bridgeType, bridgeSeq, bridgeLen, family, srcId = self.supports_unaligned_bridge(index, viralDb, outDir)

        ## 2. If bridge not found search for supplementary alignment supporting a bridge (algorithm 2)
        # If there is viralDb, dont look for aligned bridges (since viruses are not supposed to be aligned in the genome)
        if self.bridge is False and not viralDb:
            self.bridge, supportType, bridgeType, bridgeSeq, bridgeLen, family, srcId = self.supports_aligned_bridge(maxBridgeLen, minMatchPerc, minReads, minPercReads, annotations)

        self.bridgeInfo['supportType'] = supportType
        self.bridgeInfo['bridgeType'] = bridgeType
        self.bridgeInfo['bridgeSeq'] = bridgeSeq
        self.bridgeInfo['bridgeLen'] = bridgeLen
        self.bridgeInfo['family'] = family
        self.bridgeInfo['srcId'] = srcId

    def supports_unaligned_bridge(self, index, viralDb, outDir):
        '''
        Assess of supplementary cluster supports an unaligned bridge sequence. 
        In this case there will be an unaligned piece of sequence between read level
        breakpoints for clipping and supplementary alignmnet BND junction   

        Input:
            1. index: minimap2 index for consensus retrotransposon sequences + source element downstream regions
            2. outDir: output directory     

        Output: 
            1. bridge: boolean specifying if unaligned bridge found (True) or not (False)
            2. supportType: 'unaligned' or None
            3. bridgeType: type of bridge (solo, partnered or orphan)
            4. bridgeSeq: bridge sequence
            5. bridgeLen: bridge length
            6. family: retrotransposon family
            7. srcId: source element identifier
        '''    
        bridge = False 
        supportType = None
        bridgeType = None  
        bridgeSeq = None
        bridgeLen = None
        family = None 
        srcId = None 

        ## 1. Remove None events
        insertSizes = [event.insertSize for event in self.events if event.insertSize] 

        ## 2. Supported by less than X events -> Abort, no bridge
        nbEvents = len(insertSizes)

        if nbEvents < 2:
            return bridge, supportType, bridgeType, bridgeSeq, bridgeLen, family, srcId

        ## 3. Assess bridge size consistency between supporting events
        # AND select representative event if bridge found
        meanLen = np.mean(insertSizes)
        std = np.std(insertSizes)
        cv = std / meanLen * 100 

        ## Consistent bridge sizes -> call bridge 
        if cv <= 40:

            ## Compute median bridge size
            medianSize = np.median(insertSizes)

            ## Select event closest to the median size as representative
            closestSize = min(insertSizes, key=lambda x:abs(x-medianSize))

            for event in self.events:

                if event.insertSize == closestSize:
                    self.representative = event
                    break

            ## 4. Infer bridge structure 
            ## Generate fasta containing representative inserted sequence at BND junction
            fasta = formats.FASTA()
            fasta.seqDict['insert'] = self.representative.insertSeq
            insertPath = outDir + '/insert.fa'
            fasta.write(insertPath)

            ## Infer structure
            # Look for retrotransposons if there is not viralDb.
            if not viralDb:
                structure = retrotransposons.retrotransposon_structure(insertPath, index, outDir)

                # a) Resolved structure
                if ('INS_TYPE' in structure) and (structure['INS_TYPE'] is not 'unknown') and ('PERC_RESOLVED' in structure) and (structure['PERC_RESOLVED'] >= 60):
                    bridge = True
                    supportType = 'unaligned'
                    bridgeType = structure['INS_TYPE']  
                    bridgeSeq = self.representative.insertSeq
                    bridgeLen = self.representative.insertSize
                    family = ','.join(structure['FAMILY']) 
                    srcId = ','.join(structure['CYTOBAND']) if ('CYTOBAND' in structure and structure['CYTOBAND']) else None

                # b) Unresolved structure
                else: 
                    bridge = True
                    supportType = 'unaligned'
                    bridgeType = 'unknown'
                    bridgeSeq = self.representative.insertSeq
                    bridgeLen = self.representative.insertSize
                    family = None
                    srcId = None

            # Look for viruses if there is viralDb.
            else:
                structure = virus.virus_structure(insertPath, viralDb, outDir)


                # a) Resolved structure
                # TEMP
                if ('INS_TYPE' in structure) and (structure['INS_TYPE'] is not 'unknown') and ('PERC_RESOLVED' in structure) and (structure['PERC_RESOLVED'] >= 0):
                #if ('INS_TYPE' in structure) and (structure['INS_TYPE'] is not 'unknown') and ('PERC_RESOLVED' in structure) and (structure['PERC_RESOLVED'] >= 60):
                    bridge = True
                    supportType = 'unaligned'
                    bridgeType = structure['INS_TYPE']
                    bridgeSeq = self.representative.insertSeq
                    bridgeLen = self.representative.insertSize
                    family = ','.join(structure['FAMILY']) 
                    ## TODO: add virusDsc to output
                    srcId = ','.join(structure['CYTOBAND']) if ('CYTOBAND' in structure and structure['CYTOBAND']) else None

                # b) Unresolved structure
                else: 
                    bridge = True
                    supportType = 'unaligned'
                    bridgeType = 'unknown'
                    bridgeSeq = self.representative.insertSeq
                    bridgeLen = self.representative.insertSize
                    family = None
                    srcId = None

        return bridge, supportType, bridgeType, bridgeSeq, bridgeLen, family, srcId    

    def supports_aligned_bridge(self, maxBridgeLen, minMatchPerc, minReads, minPercReads, annotations):
        '''
        Assess of supplementary cluster supports a bridge sequence aligning 
        over an annotated L1 or transduced region on the reference genome

        Input:
            1. maxBridgeLen: maximum supplementary cluster length to search for a bridge
            2. minMatchPerc: minimum percentage of the supplementary cluster interval to match in a transduction or repeats database to make a bridge call
            3. minReads: minimum number of reads supporting the bridge
            4. minPercReads: minimum percentage of clipping cluster supporting reads composing the bridge
            5. annotations: dictionary containing one key per type of annotation loaded and bin databases containing annotated features as values (None for those annotations not loaded)

        Output: 
            1. bridge: boolean specifying if unaligned bridge found (True) or not (False)
            2. supportType: 'aligned' or None
            3. bridgeType: type of bridge (solo, partnered or orphan)
            4. bridgeLen: bridge length
            5. family: retrotransposon family
            6. srcId: source element identifier        
        '''
        bridge = False 
        supportType = None
        bridgeType = None  
        bridgeSeq = None        
        bridgeLen = None
        family = None 
        srcId = None 

        # Cluster interval length within limit
        # NOTE: request clusterIndex == 0
        if self.length() <= 10000:

            ## 1. Assess if supplementary cluster spans a transduced region
            if ('TRANSDUCTIONS' in annotations) and (self.ref in annotations['TRANSDUCTIONS']):
                        
                ## Select transduction bin database for the corresponding ref            
                tdBinDb = annotations['TRANSDUCTIONS'][self.ref]

                ## Intersect supplementary cluster interval with transducted regions  
                tdMatches = tdBinDb.collect_interval(self.beg, self.end, 'ALL')

                ## Sort matches in decreasing match % order
                sortedTdMatches = sorted(tdMatches, key=itemgetter(2), reverse=True)
            
                ## Check if cluster matches a transduced area (use match with longest %)
                if sortedTdMatches and sortedTdMatches[0][2] > minMatchPerc:
                    tdMatch = sortedTdMatches[0]
                
                else:
                    tdMatch = None

                ## Stop if suppl. cluster matches transduced area
                if tdMatch is not None:

                    bridge = True
                    supportType = 'aligned'
                    bridgeType = 'transduced'
                    bridgeLen = self.length()
                    family = 'L1'
                    srcId = tdMatch[0].optional['cytobandId']

            ## 2. Assess if supplementary cluster spans an annotated repeat
            if ('REPEATS' in annotations) and (self.ref in annotations['REPEATS']):

                ## Select repeats bin database for the corresponding ref
                repeatsBinDb = annotations['REPEATS'][self.ref]

                ## Intersect supplementary cluster interval with annotated repeats 
                repeatMatches = repeatsBinDb.collect_interval(self.beg, self.end, 'ALL')

                ## Sort matches in decreasing match % order
                sortedRepeatMatches = sorted(repeatMatches, key=itemgetter(2), reverse=True)

                ## Check if cluster matches a repeat (use match with longest %)
                if sortedRepeatMatches and sortedRepeatMatches[0][2] > minMatchPerc:
                    repeatMatch = sortedRepeatMatches[0]
                else:
                    repeatMatch = None

                ## Cluster matches annotated repeat area
                if repeatMatch is not None:
                    bridge = True
                    supportType = 'aligned'
                    bridgeType = 'solo'
                    bridgeLen = self.length()
                    family = repeatMatch[0].optional['family']

        return bridge, supportType, bridgeType, bridgeSeq, bridgeLen, family, srcId

    def is_junction_partner(self, bridgeReadIds, minPercReads):
        '''
        Assess if supplementary cluster is a BND junction partner

        Input:
            1. bridgeReadIds: list of read ids supporting bridge at the BND junction
            2. minPercReads: minimum perc of cluster supporting reads also supporting the junction
        
        Output:
            1. partner: boolean (True: junction partner; False: not partner)
        '''
        ## 1. Count the number of suppl. events supporting a BND junction partnership
        nbPartnerSuppl = 0

        # For each supplementary event composing the cluster 
        for supplementary in self.events:

            ## a) Read supporting also a bridge
            if supplementary.readName in bridgeReadIds:

                if supplementary.readIndex == 1:
                    nbPartnerSuppl +=1

            ## b) Read not supporting a bridge
            else:
                if supplementary.readIndex == 0:
                    nbPartnerSuppl +=1
                
        ## 2. Assess if enough % of events supporting the junction
        nbSuppl = float(len(self.events))
        percPartnered = nbPartnerSuppl / nbSuppl * 100

        # a) Make bridge call
        if percPartnered >= minPercReads:
            partner = True

        # b) No bridge call
        else:
            partner = False
        
        return partner

class DISCORDANT_cluster(cluster):
    '''
    Discordant cluster subclass
    '''
    def __init__(self, events):

        cluster.__init__(self, events, 'DISCORDANT')

        self.matesCluster = None
        self.identity = None
        self.element = self.events[0].element

        if all (event.orientation == 'PLUS' for event in events):
            self.orientation = 'PLUS'
        elif all (event.orientation == 'MINUS' for event in events):
            self.orientation = 'MINUS'
        else:
            self.orientation = 'RECIPROCAL'

        # cluster.__init__(self, events, 'DISCORDANT')

        # self.matesCluster = None
        # self.identity = self.events[0].identity
        # self.orientation = self.events[0].orientation
        # self.element = self.events[0].element
            
        # if all (event.orientation == 'PLUS' for event in events):
        #     self.orientation = 'PLUS'
        # elif all (event.orientation == 'MINUS' for event in events):
        #     self.orientation = 'MINUS'
        # else:
        #     self.orientation = 'RECIPROCAL'
    
    def setIdentity(self):
        '''
        Set identity based on events identity
        '''
        identities = set([event.identity for event in self.events])
           
        if len(identities) == 1:
            self.identity = self.events[0].identity
        else:
            self.identity = list(identities)
            
    def mates_start_interval(self):
        '''
        Compute mates start alignment position interval. 
        
        Only makes senses if all the mates composing the cluster are clustered 
        in the same chromosome

        Output:
            1. beg: interval begin position
            2. end: interval end position
        '''
        starts = [event.mateStart for event in self.events]
        beg = min(starts)
        end = max(starts)
        return beg, end
        
    def create_matesCluster(self):
        '''
        Create discordant read pair cluster for mates
        '''
        ## 1. Produce discordant objects for mates:
        mates = events.discordants2mates(self.events)

        ## 2. Create discordant cluster for mates
        matesCluster = DISCORDANT_cluster(mates)
        
        return matesCluster

class META_cluster():
    '''
    Meta cluster class
    '''
    number = 0 # Number of instances

    def __init__(self, clusters):
        '''
        '''
        META_cluster.number += 1 # Update instances counter
        self.id = 'META_' + str(META_cluster.number)

        # Define list of events composing the cluster 
        self.events = list(itertools.chain(*[cluster.events for cluster in clusters]))

        # Set cluster's reference, begin and end position
        self.ref, self.beg, self.end = self.coordinates()
        self.refLeftBkp = None
        self.refRightBkp = None
        
        # Organize events into subclusters
        self.subclusters = self.create_subclusters()
        # NOTE MERGE SR2020: To avoid key META error in clustering.reciprocal_overlap_clustering
        self.rawSubclusters = clusters

        # Set some metacluster properties as None
        self.bkpPos = None
        self.mutOrigin = None
        self.failedFilters = None
        self.consensusEvent = None                
        self.insertHits = None
        self.nbTotal, self.nbTumour, self.nbNormal, self.nbINS, self.nbDEL, self.nbCLIPPING = [None, None, None, None, None, None] 
        self.nbReadsTotal, self.nbReadsTumour, self.nbReadsNormal, self.reads, self.readsTumour, self.readsNormal = [None, None, None, None, None, None] 
        self.cv = None
        self.repreLeftSeq = None
        self.repreRightSeq = None
        self.consLeftSeq = None
        self.consRightSeq = None
        self.intLeftBkp = None
        self.intRightBkp = None
        self.rightClipType = None
        self.leftClipType = None

        ## Short reads:
        self.orientation = None
        # ins type and identity
        self.ins_type = None
        self.identity = None
        self.plus_id = None
        self.minus_id = None
        # transduction attributes
        self.src_id = None
        self.plus_src_id = None
        self.minus_src_id = None
        self.src_end = None
        self.plus_src_end = None
        self.minus_src_end = None
        self.src_type = None
        # polyA
        self.pA = None
        self.plus_pA = None
        self.minus_pA = None
        # mei characteristics
        self.TSD = None
        self.strand = None
        self.repeatAnnot = None

        # Short reads:
        if hasattr(self.events[0], 'identity'):
            self.identity = self.events[0].identity
            
        if hasattr(clusters[0], 'orientation'):
            if all (cluster.orientation == 'PLUS' for cluster in clusters):
                self.orientation = 'PLUS'
            elif all (cluster.orientation == 'MINUS' for cluster in clusters):
                self.orientation = 'MINUS'
            else:
                self.orientation = 'RECIPROCAL'
                    
        # Update input cluster's clusterId attribute
        for cluster in clusters:
            cluster.clusterId = self.id

        # Initialize dictionaries 
        self.SV_features = {}
        self.supplClusters = {}
        self.bridge = None
    
    def identity_shortReads_ME(self):
        '''
        Define plus and minus cluster identities, and set polyA attributes if polyA present
        '''
        self.identity = None
        self.plus_id = None
        self.minus_id = None

        ## 1. collect events identities
        plus_identities = [event.identity for event in self.events if event.orientation == 'PLUS']
        minus_identities = [event.identity for event in self.events if event.orientation == 'MINUS']     
        
        # get total counts
        countPlus = len(plus_identities)
        countMinus = len(minus_identities)
        
        ## 2. set pA support
        if 'Simple_repeat' in plus_identities:
            self.plus_pA, self.pA = True, True
            plus_identities = [value for value in plus_identities if value != 'Simple_repeat']
        
        if 'Simple_repeat' in minus_identities:
            self.minus_pA, self.pA = True, True
            minus_identities = [value for value in minus_identities if value != 'Simple_repeat']
                    
        ## 3. set clusters identity
        plus_identityDict = Counter(plus_identities)
        minus_identityDict = Counter(minus_identities)
        
        # select most supported identity for plus cluster
        if plus_identityDict:
            
            # select identity with the max support
            plus_identityDict_max = max(plus_identityDict.items(), key=itemgetter(1))
            
            # if the support is > 10% of total identity counts
            if plus_identityDict_max[1]/countPlus > 0.1:
                
                plus_id = plus_identityDict_max[0].split('_')
                self.plus_id = plus_id[0]
                
                # if it's a TD
                if len(plus_id) == 3:
                    self.plus_src_id, self.plus_src_end = plus_id[1], plus_id[2]
            
        # if polyA ('Simple_repeat') support is more than 90%
        if not self.plus_id and countPlus > 0:
            self.plus_id = 'Simple_repeat'
        
        # select most supported identity for minus cluster    
        if minus_identityDict:
            
            # select identity with the max support
            minus_identityDict_max = max(minus_identityDict.items(), key=itemgetter(1))
            
            # if the support is > 10% of total identity counts
            if minus_identityDict_max[1]/countMinus > 0.1:
                
                minus_id = minus_identityDict_max[0].split('_')
                self.minus_id = minus_id[0]
                
                # if it's a TD
                if len(minus_id) == 3:
                    self.minus_src_id, self.minus_src_end = minus_id[1], minus_id[2]
            
        # if polyA ('Simple_repeat') support is more than 90%
        if not self.minus_id and countMinus > 0:
            self.minus_id = 'Simple_repeat'
                   
    def sort(self):
        '''
        Sort events in increasing coordinates order
        '''
        self.events.sort(key=lambda event: event.beg)

    def coordinates(self):
        '''
        Compute cluster ref, beg and end coordinates. 
        
        Begin and end will correspond to the left and rightmost positions, respectively
        '''
        # Sort events from lower to upper beg coordinates
        self.sort()  

        # Define cluster coordinates 
        ref = self.events[0].ref
        beg = min([event.beg for event in self.events])
        end = max([event.end for event in self.events])
        
        return ref, beg, end

    def mean_pos(self):
        '''
        Compute cluster mean genomic position and confidence interval around the mean 
        
        Output:
            1. pos: mean genomic position
            2. cipos: confidence interval around the mean position
        '''
        ## 1. Collect all begin positions
        begs = [event.beg for event in self.events]

        ## 2. Compute mean position and standard error 
        mean = int(np.mean(begs))
        sem = scipy.stats.sem(begs)

        ## 3. Compute confidence interval positions
        if (sem == 0) or (math.isnan(sem)):
            CI = (mean, mean)
        else:
            CI = scipy.stats.t.interval(0.95, len(begs)-1, loc=mean, scale=sem)
            
        ## Make confidence interval relative to the mean position
        CI = (int(CI[0] - mean), int(CI[1] - mean))
        
        return mean, CI

    def create_subclusters(self):
        '''
        Organize cluster composing events into subclusters

        Output:
            1. subclusters: dictionary containing cluster types as keys and cluster object as values
        '''
        ## 1. Separate events according to their type into multiple lists ##
        eventTypes = events.separate(self.events)

        ## 2. Create subclusters ##
        subclusters = {}

        for eventType, eventList in eventTypes.items():

            ## Create subcluster
            subcluster = create_cluster(eventList, eventType) 

            ## Set subcluster metacluster id attribute:
            subcluster.clusterId = self.id

            ## Add subcluster to the dict
            subclusters[eventType] = subcluster 

        return subclusters

    def add(self, clusters2add):
        '''
        Add a list of clusters to the metacluster 

        Input:
            1. clusters2add: List of clusters to be added 
        '''
        ## 0. Update metacluster's id attribute
        for cluster in clusters2add:
            cluster.clusterId = self.id

        ## 1. Add events within the input cluster to the metacluster ##
        events2add = list(itertools.chain(*[cluster.events for cluster in clusters2add]))
        self.events = self.events + events2add

        ## 2. Resort and redefine metacluster begin and end coordinates ##
        self.ref, self.beg, self.end = self.coordinates()

        ## 3. Separate events according to their type into multiple lists ##
        eventTypes = events.separate(events2add)

        ## 4. Add events to the subclusters ##
        for eventType, eventList in eventTypes.items():
            
            # a) Create subcluster if not pre-existing one
            if eventType not in self.subclusters:
         
                ## Create subcluster
                subcluster = create_cluster(eventList, eventType) 
            
                ## Set subcluster metacluster id attribute:
                subcluster.clusterId = self.id

                ## Add subcluster to the dict
                self.subclusters[eventType] = subcluster 

            # b) Add events to pre-existing subcluster
            else:
                self.subclusters[eventType].add(eventList)

        # Update input cluster's clusterId attribute
        for cluster in clusters2add:
            cluster.clusterId = self.id
        # NOTE MERGE SR2020: To avoid key META error in clustering.reciprocal_overlap_clustering
        # Also add clusters to rawSubclusters:
        self.rawSubclusters.extend(clusters2add)

    def addEvents(self, eventsList):
        '''

        Input:
            1. events: List of events to be added to the metacluster
        '''
        ## 1. Add events to the cluster ##
        previous = self.events
        self.events = self.events + eventsList

        ## 2. Resort and redefine cluster begin and end coordinates ##
        self.ref, self.beg, self.end = self.coordinates()

        ## 3. Separate events according to their type into multiple lists ##
        eventTypes = events.separate(eventsList)

        ## 4. Add events to the subclusters ##
        for eventType, eventList in eventTypes.items():
            
            # a) Create subcluster if not pre-existing one
            if eventType not in self.subclusters:
         
                ## Create subcluster
                subcluster = create_cluster(eventList, eventType) 
            
                ## Add subcluster to the dict
                self.subclusters[eventType] = subcluster 

            # b) Add events to pre-existing subcluster
            else:
                self.subclusters[eventType].add(eventList)

    def remove(self, events2remove):
        '''
        Remove a list of events from the metacluster and corresponding subclusters

        Input:
            1. events2remove: List of events to be removed 
        '''
        ## 1. Remove events from the metacluster ##
        self.events = [event for event in self.events if event not in events2remove]

        ## 2. Resort and redefine metacluster begin and end coordinates ##
        if self.events:
            self.ref, self.beg, self.end = self.coordinates()

        ## 3. Separate events according to their type into multiple lists ##
        eventTypes = events.separate(events2remove)

        ## 4. Remove events from the subclusters ##
        for eventType, eventList in eventTypes.items():

            if eventType in self.subclusters:
                self.subclusters[eventType].remove(eventList)
                
            else:
                log.info('WARNING at remove method from META_cluster class. Event with unkown type')

    def collect_reads(self):
        '''
        Create FASTA object containing metacluster supporting reads.
        
        Output:
            1. FASTA: FASTA object containing metacluster supporting reads
        '''
        ## Initiate FASTA object
        FASTA = formats.FASTA()

        ## For each event composing the metacluster         
        for event in self.events: 

            ## A) Read sequence supporting the event not included in the FASTA yet
            if event.readName not in FASTA.seqDict:
                FASTA.seqDict[event.readName] = event.readSeq 

            ## B) Read sequence supporting the event already included in the FASTA
            else:

                # Current read sequence aligned as primary -> replace previously included sequence (therefore supplementary)
                if not event.supplementary:
                    FASTA.seqDict[event.readName] = event.readSeq 
     
        return FASTA

    def supportingReads(self):
        '''
        Compute the total number of metacluster supporting reads and generate a list of supporting read ids

        Output:
            1. nbTotal: Total number of metacluster supporting reads
            2. nbTumour: Number of metacluster supporting reads in the tumour 
            3. nbNormal: Number of metacluster supporting reads in the normal
            4. reads: List containing all the metacluster supporting reads
            5. readsTumour: List containing metacluster supporting reads in the tumour sample
            6. readsNormal: List containing metacluster supporting reads in the normal sample
        '''
        reads = []
        readsTumour = []
        readsNormal = []

        ## 1. Create non-redundant list of metacluster supporting reads
        for event in self.events: 

            # Add read name not already included in the list  
            if event.readName not in reads:
                reads.append(event.readName)

            ## Tumour and matched normal
            # a) Event identified in the TUMOUR sample
            if event.sample == "TUMOUR":
            
                # Add read name not already included in the list  
                if event.readName not in readsTumour:
                    readsTumour.append(event.readName)
            
            # b) Event identified in the matched NORMAL sample
            elif event.sample == "NORMAL":

                # Add read name not already included in the list  
                if event.readName not in readsNormal:
                    readsNormal.append(event.readName)            

        ## 2. Compute number of supporting reads
        nbTotal = len(reads) 

        # a) Unpaired mode
        if self.events[0].sample is None:
            nbTumour, nbNormal, readsTumour, readsNormal = [None, None, None, None]

        # b) Paired mode
        else:

            nbTumour = len(readsTumour) # tumour
            nbNormal = len(readsNormal) # normal

        return nbTotal, nbTumour, nbNormal, reads, readsTumour, readsNormal

    def nbEvents(self):
        '''
        Return the number of events composing the metacluster. 
        '''
        ## Initialize counters
        nbTumour = 0
        nbNormal = 0

        nbINS = 0
        nbDEL = 0
        nbCLIPPING = 0   

        # For each event composing the metacluster
        for event in self.events:

            ## Tumour and matched normal counts
            # a) Event identified in the TUMOUR sample
            if event.sample == "TUMOUR":
                nbTumour += 1
            
            # b) Event identified in the matched NORMAL sample
            elif event.sample == "NORMAL":
                nbNormal += 1
            
            # c) SINGLE sample mode
            else:
                nbTumour = None
                nbNormal = None
                            
            ## Event type counts
            # a) INS event
            if event.type == 'INS':
                nbINS += 1

            # b) DEL event
            elif event.type == 'DEL':
                nbDEL += 1

            # c) CLIPPING event
            else:
                nbCLIPPING += 1            

        nbTotal = len(self.events)

        return nbTotal, nbTumour, nbNormal, nbINS, nbDEL, nbCLIPPING

    def percDuplicates(self):
        '''
        Return the number of events that compose the metacluster and are labbeled as duplicate in the input bam file. 
        '''
        ## Initialize counters
        nbDiscordantDuplicatesTumour = 0
        nbDiscordantDuplicatesNormal = 0
        nbDiscordantDuplicatesTotal = 0
        nbClippingDuplicatesTumour = 0
        nbClippingDuplicatesNormal = 0
        nbClippingDuplicatesTotal = 0

        nbDiscordantTumour = 0
        nbClippingTumour = 0
        nbDiscordantNormal = 0
        nbClippingNormal = 0
        nbDiscordantTotal = 0
        nbClippingTotal = 0

        # For each event composing the metacluster
        for event in self.events:

            ## Tumour and matched normal counts
            # a) Event identified in the TUMOUR sample
            if event.sample == "TUMOUR":
                if event.isDup == True:
                    if event.type == 'DISCORDANT':
                        nbDiscordantDuplicatesTumour += 1
                    elif event.type == 'CLIPPING':
                        nbClippingDuplicatesTumour  += 1
                elif event.isDup == False:
                    if event.type == 'DISCORDANT':
                        nbDiscordantTumour += 1
                    elif event.type == 'CLIPPING':
                        nbClippingTumour  += 1

            elif event.sample == "NORMAL":
                if event.isDup == True:
                    if event.type == 'DISCORDANT':
                        nbDiscordantDuplicatesNormal += 1
                    elif event.type == 'CLIPPING':
                        nbClippingDuplicatesNormal  += 1
                elif event.isDup == False:
                    if event.type == 'DISCORDANT':
                        nbDiscordantNormal += 1
                    elif event.type == 'CLIPPING':
                        nbClippingNormal  += 1
            
            # c) SINGLE sample mode
            elif event.sample == None:
                if event.isDup == True:
                    if event.type == 'DISCORDANT':
                        nbDiscordantDuplicatesTotal += 1
                    elif event.type == 'CLIPPING':
                        nbClippingDuplicatesTotal += 1
                elif event.isDup == False:
                    if event.type == 'DISCORDANT':
                        nbDiscordantTotal += 1
                    elif event.type == 'CLIPPING':
                        nbClippingTotal  += 1
                            
        if nbDiscordantDuplicatesTotal == 0:
            nbDiscordantDuplicatesTotal = nbDiscordantDuplicatesTumour + nbDiscordantDuplicatesNormal
            nbClippingDuplicatesTotal = nbClippingDuplicatesTumour + nbClippingDuplicatesNormal
            nbDiscordantTotal = nbDiscordantTumour + nbDiscordantNormal
            nbClippingTotal = nbClippingTumour + nbClippingNormal

        if nbDiscordantDuplicatesTotal > 0:
            percDiscordantDuplicates = (nbDiscordantDuplicatesTotal / (nbDiscordantDuplicatesTotal + nbDiscordantTotal)) * 100
        else:
            percDiscordantDuplicates = 0
        if nbClippingDuplicatesTotal > 0:
            percClippingDuplicates = (nbClippingDuplicatesTotal / (nbClippingDuplicatesTotal + nbClippingTotal)) * 100
        else:
            percClippingDuplicates = 0


        return percDiscordantDuplicates, percClippingDuplicates

    def nbDISCORDANT(self):
        '''
        Return the number of discordant events composing the metacluster. 
        '''                       
        nbDISCORDANT = len(set([event.readName for event in self.events if event.type == 'DISCORDANT']))

        return nbDISCORDANT

    def nbSUPPLEMENTARY(self):
        '''
        Return the number of supplementary events composing the metacluster. 
        '''        
        nbSUPPL = len(set([event.readName for event in self.events if event.type == 'SUPPLEMENTARY']))
        
        return nbSUPPL

    def nbCLIPPINGS(self):
        '''
        Return the number of clipping events composing the metacluster. 
        '''        
        nbCLIPP = len(set([event.readName for event in self.events if event.type == 'CLIPPING']))
        
        return nbCLIPP

    def supportingCLIPPING(self, buffer, confDict, bam, normalBam, mode):
        # Note: This function works but you have to allow duplicates in the clipping 

        # Make custom conf. dict for only selecting duplicates
        clippingConfDict = dict(confDict)
        clippingConfDict['targetSV'] = ['CLIPPING']
        clippingConfDict['minMAPQ'] = 0

        clippingEventsDict = {}

        ## Define region
        binBeg = self.beg - buffer if self.beg > buffer else 0
        
        # TODO check as above
        binEnd = self.end

        ref = self.ref

        if mode == "SINGLE":
            clippingEventsDict = bamtools.collectSV(ref, binBeg, binEnd, bam, clippingConfDict, None)
        elif mode == "PAIRED":
            clippingEventsDict = bamtools.collectSV_paired(ref, binBeg, binEnd, bam, normalBam, clippingConfDict)

        ## When the discordant cluster is RIGHT, add the biggest right clipping cluster if any:
        if all (event.side == 'PLUS' for event in self.events):
            
            ## Get clipping clusters:
            clippingRightEventsDict = dict((key,value) for key, value in clippingEventsDict.items() if key == 'RIGHT-CLIPPING')
            CLIPPING_cluster = self.add_clippingEvents(ref, binBeg, binEnd, clippingRightEventsDict, ['RIGHT-CLIPPING'], confDict)

        ## When the discordant cluster is LEFT, add the biggest left clipping cluster if any:
        elif all (event.side == 'MINUS' for event in self.events):
            
            ## Get clipping clusters:
            clippingLeftEventsDict = dict((key,value) for key, value in clippingEventsDict.items() if key == 'LEFT-CLIPPING')
            CLIPPING_cluster = self.add_clippingEvents(ref, binBeg, binEnd, clippingLeftEventsDict, ['LEFT-CLIPPING'], confDict)

        # TODO if it is reciprocal
        else:
            CLIPPING_cluster = self.add_clippingEvents(ref, binBeg, binEnd, clippingEventsDict, ['RIGHT-CLIPPING', 'LEFT-CLIPPING'], confDict)

        return CLIPPING_cluster

    def supportingCLIPPING_sonia(self, confDict, bam, normalBam):
        '''
        Collect supporting clippings 
        '''
        
        # Make custom conf. dict
        clippingConfDict = dict(confDict)
        clippingConfDict['targetEvents'] = ['CLIPPING']
        clippingConfDict['minMAPQ'] = 1
        clippingConfDict['readFilters'] = ['SMS']
        clippingConfDict['minClusterSize'] = 1

        clippingEventsDict = {}

        # Define region
        ref = self.ref
        beg = self.refLeftBkp if self.refLeftBkp is not None else self.beg
        end = self.refRightBkp if self.refRightBkp is not None else self.end

        # Collect clippings
        if not normalBam:
            clippingEventsDict = bamtools.collectSV(ref, beg, end, bam, clippingConfDict, None)
            
        else:
            clippingEventsDict = bamtools.collectSV_paired(ref, beg, end, bam, normalBam, clippingConfDict)
        
        # When cluster orientation is 'PLUS', add the biggest right clipping cluster if any:
        if self.orientation == 'PLUS':
            
            ## Get clipping clusters:
            clippingRightEventsDict = dict((key,value) for key, value in clippingEventsDict.items() if key == 'RIGHT-CLIPPING')
            CLIPPING_cluster = self.add_clippingEvents(ref, beg, end, clippingRightEventsDict, ['RIGHT-CLIPPING'], confDict)

        # When cluster orientation is LEFT, add the biggest left clipping cluster if any:
        elif self.orientation == 'MINUS':
            
            ## Get clipping clusters:
            clippingLeftEventsDict = dict((key,value) for key, value in clippingEventsDict.items() if key == 'LEFT-CLIPPING')
            CLIPPING_cluster = self.add_clippingEvents(ref, beg, end, clippingLeftEventsDict, ['LEFT-CLIPPING'], confDict)

        # If it is reciprocal
        elif self.orientation == 'RECIPROCAL':
            
            CLIPPING_cluster = self.add_clippingEvents(ref, beg, end, clippingEventsDict, ['RIGHT-CLIPPING', 'LEFT-CLIPPING'], confDict)

        return CLIPPING_cluster
        
    def add_clippingEvents(self, ref, binBeg, binEnd, clippingEventsDict, eventTypes, confDict):
        '''
        Create clipping clusters and add them to metacluster
        '''
        binSizes = [1000]
        clippingBinDb = structures.create_bin_database_interval(ref, binBeg, binEnd, clippingEventsDict, binSizes)
        binSize = clippingBinDb.binSizes[0]
        CLIPPING_clusters = clustering.distance_clustering(clippingBinDb, binSize, eventTypes, 'CLIPPING', confDict['maxBkpDist'], confDict['minClusterSize']) 

        print('clippingBinDb')
        print(clippingBinDb)
        
        # If there is a clipping cluster
        for CLIPPING_cluster in CLIPPING_clusters:
            
            ## Add cluster's reads to the discordant metacluster:
            self.addEvents(CLIPPING_cluster.events)

        return CLIPPING_clusters
 
    def polish(self, confDict, reference, outDir):
        '''
        Polish metacluster consensus sequence

        Input: 
            1. confDict: 
                * technology     -> sequencing technology (NANOPORE, PACBIO or ILLUMINA)
                * rounds         -> number of polishing rounds to be attempled. 0 means no polishing
            
            2. reference: path to reference genome in fasta format    
            3. outDir: Output directory
        
        Output: Update 'consensusFasta' attribute with the polished sequence
        '''
        ## 0. Create directory 
        unix.mkdir(outDir)        

        ## 1. Use cluster supporting reads to polish metacluster consensus sequence ##
        ## 1.1 Write raw consensus sequence
        unpolishedFasta = outDir + '/raw_consensus.fa'
        self.consensusFasta.write(unpolishedFasta)

        ## 1.2 Collect reads that will be used to polish the consensus sequence        
        supportingReads = self.collect_reads()

        ## Remove consensus from FASTA
        consensusReadName = list(self.consensusFasta.seqDict.keys())[0]

        if consensusReadName in supportingReads.seqDict:
            del supportingReads.seqDict[consensusReadName]

        ## Write supporting reads FASTA 
        supportingReadsFasta = outDir + '/supportingReads.fa'
        supportingReads.write(supportingReadsFasta)
            
        ## 2. Consensus polishing 
        polishedFasta = assembly.polish_racon(unpolishedFasta, supportingReadsFasta, confDict['technology'], confDict['rounds'], outDir)

        ## If polishing is successful replace consensus by new polished sequence
        if polishedFasta is not None:
            polished = formats.FASTA()
            polished.read(polishedFasta)     
            self.consensusFasta = polished

        return 

    def consensus_event(self, confDict, reference, offset, outDir):
        '''
        Define metacluster´s consensus event based on consensus sequence realignment

        Input: 
            1. confDict: 
                * targetSV       -> list with target SV (INS: insertion; DEL: deletion; CLIPPING: left and right clippings)
                * minMAPQ        -> minimum mapping quality
                * minCLIPPINGlen -> minimum clipping lenght
                * minINDELlen    -> minimum INS and DEL lenght
                * overhang       -> Number of flanking base pairs around the INDEL events to be collected from the supporting read. If 'None' the complete read sequence will be collected)            2. outDir: output directory
            
            2. reference: path to reference genome in fasta format    
            3. offset: number of base pairs to extend cluster begin and end coordinates when defining target region for consensus sequence realignment
            4. outDir: Output directory
        
        Output: Update 'consensusEvent' attribute in the metacluster class
        '''

        ## 1. Local realignment of the consensus sequence into the SV genomic interval ##
        ## 1.1 Write consensus sequence into a fasta file 
        consensusFile = outDir + '/consensus.fa'
        self.consensusFasta.write(consensusFile)

        ## 1.2 Define SV cluster genomic interval
        # ------------------<***SV_cluster***>-----------------
        #       <--offset-->                  <--offset-->
        intervalBeg = self.beg - offset
        intervalBeg = intervalBeg if intervalBeg >= 0 else 0 ## Set lower bound
        intervalEnd = self.end + offset
        intervalCoord = self.ref + ':' + str(intervalBeg) + '-' + str(intervalEnd)
            
        ## 1.3 Do realignment
        # ------------------<***SV_cluster***>-----------------
        #         -------------consensus_seq-------------
        BAM = alignment.targeted_alignment_minimap2(consensusFile, intervalCoord, reference, outDir, 'BAM')
 
        ## Continue if realignment is succesfull 
        if BAM is not None:
            
            ## 2. Search for metacluster´s consensus event  ##
            ## 2.1 Extract events from consensus sequence realignment
            # ------------------<***SV_cluster***>-----------------
            #        -------------consensus_seq-------------
            #               <--->----------------<---> overhang (100 bp)
            #                   event_search_space         
            overhang = 100
            clusterIntervalLen = self.end - self.beg
            targetBeg = offset - overhang
            targetEnd = offset + clusterIntervalLen + overhang            
            eventsDict = bamtools.collectSV(intervalCoord, targetBeg, targetEnd, BAM, confDict, None)

            ## 2.2 Define consensus event based on the events resulting from consensus sequence realignment
            ## A) Metacluster supports an INS and realignment leads to one INS event 
            if (self.SV_type == 'INS') and ('INS' in eventsDict) and (len(eventsDict['INS']) == 1):

                ## Convert coordinates
                self.consensusEvent = alignment.targetered2genomic_coord(eventsDict['INS'][0], self.ref, intervalBeg)

            ## B) Metacluster supports an INS and realignment leads to multiple INS events
            elif (self.SV_type == 'INS') and ('INS' in eventsDict) and (len(eventsDict['INS']) > 1):

                ## Do merging
                merged = events.merge_INS(eventsDict['INS'])

                ## Convert coordinates
                self.consensusEvent = alignment.targetered2genomic_coord(merged, self.ref, intervalBeg)

            ## C) Metacluster supports a DEL and realignment leads to one DEL event 
            elif (self.SV_type == 'DEL') and ('DEL' in eventsDict) and (len(eventsDict['DEL']) == 1):

                ## Convert coordinates
                self.consensusEvent = alignment.targetered2genomic_coord(eventsDict['DEL'][0], self.ref, intervalBeg)

            ## D) Metacluster supports a DEL and realignment leads to multiple DEL events (TO DO)
            #elif (self.SV_type == 'DEL') and (len(eventsDict['DEL']) > 1):
                    ## Do merging
                    ## Convert coordinates
                        
            ## E) Metacluster supports an INS and realignment leads to one left and one right CLIPPING 
            elif (self.SV_type == 'INS') and ('LEFT-CLIPPING' in eventsDict) and ('RIGHT-CLIPPING' in eventsDict) and (len(eventsDict['RIGHT-CLIPPING']) == 1) and (len(eventsDict['LEFT-CLIPPING']) == 1):

                ## Compute the inserted sequence length as the difference between both clipping breakpoints in the long sequence 
                rightClipping = eventsDict['RIGHT-CLIPPING'][0]
                leftClipping = eventsDict['LEFT-CLIPPING'][0]
                length = leftClipping.readBkp - rightClipping.readBkp

                ## Create consensus INS from clipping alignments if inserted sequence found
                if length >= confDict['minINDELlen']:

                    event = events.INS(rightClipping.ref, rightClipping.beg, rightClipping.end, length, rightClipping.readName, rightClipping.readSeq, rightClipping.readBkp, None, None)
        
                    ## Convert coordinates
                    self.consensusEvent = alignment.targetered2genomic_coord(event, self.ref, intervalBeg)
    
                else:
                    print('INS_NOT_FOUND!')

            ## F) Another possibility (Don´t do anything, leave previous. Later we may need to include new conditions)

    def determine_SV_type(self, minINDELlen, technology, outDir): 
        '''
        Determine the type of structural variant (SV) supported by the metacluster and select a consensus metacluster supporting sequence and event

        SV types:
            INSERTION: 
                        >>>>>>>>>>>>>/////INS/////>>>>>>>>>>>    * Completely spanned insertion
                        >>>>>>>>>>>>>/////INS///                 * Insertion partially spanned
                                      ////INS/////>>>>>>>>>>>
            DELETION:   >>>>>>>>>>>>>-----DEL----->>>>>>>>>>>    

            DUPLICATION (TO DO) 
            INVERSION (TO DO)
            BREAK END (TO DO)      

        Input:
            1. minINDELlen: minimum INS and DEL lenght
            2. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
            3. outDir: Output directory

        Output, set the following object attributes:
            
            1. SV_type: structural variant type supported by the metacluster. None if sv type 
            2. consensusEvent: consensus event attribute 
            3. consensusFasta: fasta file containing consensus metacluster supporting sequence 
        '''
        ## Create output directory 
        unix.mkdir(outDir)
        subClusterTypes = list(self.subclusters.keys())

        ## A) Metacluster supports an insertion:
        if ('INS' in subClusterTypes) and ('DEL' not in subClusterTypes):
            self.SV_type = 'INS'

            ## Select consensus INS event and sequence
            self.consensusEvent = self.subclusters['INS'].pick_median_length()
            
            self.consensusFasta = formats.FASTA()
            self.consensusFasta.seqDict[self.consensusEvent.readName] = self.consensusEvent.readSeq

        ## B) Metacluster supports a deletion:
        elif ('DEL' in subClusterTypes) and ('INS' not in subClusterTypes):
            self.SV_type = 'DEL'

            ## Select consensus DEL event and sequence
            self.consensusEvent = self.subclusters['DEL'].pick_median_length()

            self.consensusFasta = formats.FASTA()
            self.consensusFasta.seqDict[self.consensusEvent.readName] = self.consensusEvent.readSeq
                        
        ## C) Metacluster only composed by double clipping -> Long insertion candidate
        elif (len(subClusterTypes) == 2) and ('RIGHT-CLIPPING' in subClusterTypes) and ('LEFT-CLIPPING' in subClusterTypes):

            self.consensusEvent = None                

            ## Assess if clipping clusters support an insertion
            is_INS, self.consensusFasta = double_clipping_supports_INS(self.subclusters['RIGHT-CLIPPING'], self.subclusters['LEFT-CLIPPING'], minINDELlen, technology, outDir)

            ## a) Double clippings support an INS
            if is_INS:
                self.SV_type = 'INS'

            # CHANGE 2020!
            ## b) Double clippings not support an INS, but they are still BND
            else: 
                self.SV_type = 'INS_noOverlap'

        ## D) Metacluster only composed by one clipping side (left or right) -> break end
        elif (len(subClusterTypes) == 1) and (('RIGHT-CLIPPING' in subClusterTypes) or ('LEFT-CLIPPING' in subClusterTypes)):

            self.SV_type = 'BND'

            ## Retrieve clipping cluster
            # a) Right clipping
            if 'RIGHT-CLIPPING' in subClusterTypes:
                clippingCluster = self.subclusters['RIGHT-CLIPPING']
                
            # b) Left clipping
            else:
                clippingCluster = self.subclusters['LEFT-CLIPPING']

            ## Determine BND breakpoint position
            self.bkpPos = clippingCluster.infer_breakpoint()

            ## For each clipping event composing the cluster select complementary suppl. alignments
            complementaryAlignments = clippingCluster.search4complAlignments()

            ## Cluster supplementary alignment positions
            self.supplClusters =  clippingCluster.cluster_suppl_positions(complementaryAlignments)

        ## E) Other combination -> Unknown SV type (Temporal, extend later)
        else:
            self.SV_type = None
            self.consensusEvent = None                
            self.consensusFasta = None

    def search4bridge(self, maxBridgeLen, minMatchPerc, minReads, minPercReads, annotations, index, viralDb, outDir):
        '''    
        Search for a transduction or solo repeat bridge at metacluster BND junction

        Input:
            1. maxBridgeLen: maximum supplementary cluster length to search for a bridge
            2. minMatchPerc: minimum percentage of the supplementary cluster interval to match in a transduction or repeats database to make a bridge call
            3. minReads: minimum number of reads supporting the bridge
            4. minPercReads: minimum percentage of clipping cluster supporting reads composing the bridge
            5. annotations: dictionary containing one key per type of annotation loaded and bin databases containing annotated features as values (None for those annotations not loaded)
            6. index: minimap2 index for consensus retrotransposon sequences + source element downstream regions
            7. outDir: output directory

        Output: Update 'bridgeClusters' and 'bridgeType' metacluster attributes
        '''   
        ## 1. Generate list containing all the clusters of suppl. alignments
        # Note: consider to organize the clusters directly in a list, so not conversion will be needed
        supplClusters = structures.dict2list(self.supplClusters)

        ## 2. Collect all the suppl. cluster supporting a bridge
        bridgeClusters = {}
        bridgeClusters['aligned'] = [] 
        bridgeClusters['unaligned'] = []

        ## For each cluster
        for cluster in supplClusters:

            cluster.support_bridge(maxBridgeLen, minMatchPerc, minReads, minPercReads, annotations, index, viralDb, outDir)

            ## Add cluster supporting bridge to the dictionary
            # a) Aligned bridge
            if (cluster.bridge) and (cluster.bridgeInfo['supportType'] == 'aligned'):
                bridgeClusters['aligned'].append(cluster)

            # b) Unaligned bridge
            elif (cluster.bridge) and (cluster.bridgeInfo['supportType'] == 'unaligned'):
                bridgeClusters['unaligned'].append(cluster)

        ## 3. Create bridge based on supporting suppl. alignments
        # Compute number of suppl. clusters supporting aligned and unaligned bridge
        nbAligned = len(bridgeClusters['aligned'])
        nbUnaligned = len(bridgeClusters['unaligned'])

        ## A) No bridge found if:
        # 1. No suppl. cluster supporting bridge found OR
        # 2. Multipe clusters supporting an unaligned bridge (ambigous, only one expected..)
        if (nbAligned == 0 and nbUnaligned == 0) or (nbUnaligned > 1):
            self.bridge = None

        ## B) A single cluster supporting an unaligned bridge was found 
        elif nbUnaligned == 1:

            ## Create bridge based on unaligned bridge information
            supplUnaligned = bridgeClusters['unaligned'][0]
            self.bridge = BRIDGE(bridgeClusters['unaligned'], supplUnaligned.bridgeInfo)

            ## Incorporate suppl. alignments providing alignment support to the bridge 
            # For each suppl. alignment
            for cluster in bridgeClusters['aligned']:

                ## Add suppl. alignment to bridge if it is consistent
                # a) Solo bridge
                if self.bridge.bridgeType == 'solo':

                    ## Consistent if same bridge type and family
                    if (cluster.bridgeInfo['bridgeType'] == 'solo') and (self.bridge.family == cluster.bridgeInfo['family']):
                        self.bridge.add([cluster])

                # b) Orphan bridge
                elif self.bridge.bridgeType == 'orphan':

                    ## Consistent if same bridge type and source element
                    if (cluster.bridgeInfo['bridgeType'] == 'transduced') and (self.bridge.srcId == cluster.bridgeInfo['srcId']):
                        self.bridge.add([cluster])

                # c) Partnered bridge
                elif self.bridge.bridgeType == 'partnered':

                    ## Consistent if solo and same family
                    if (cluster.bridgeInfo['bridgeType'] == 'solo') and (self.bridge.family == cluster.bridgeInfo['family']):
                        self.bridge.add([cluster])

                    ## Consistent if transduced and same source element
                    elif (cluster.bridgeInfo['bridgeType'] == 'transduced') and (self.bridge.srcId == cluster.bridgeInfo['srcId']):
                        self.bridge.add([cluster])

                # d) Unknown bridge type
                elif self.bridge.bridgeType == 'unknown':
                
                    ## Improve later to infer bridge type based on aligned info
                    self.bridge.add([cluster])

        ## C) One or multiple suppl. clusters supporting an aligned bridge
        else:
            
            ## Determine bridge type
            types = list(set([cluster.bridgeInfo['bridgeType'] for cluster in bridgeClusters['aligned']]))
            nbTypes = len(types)

            # a) Partnered
            if all(bridgeType in types for bridgeType in ['solo', 'transduced']):
                bridgeType = 'partnered'

            # b) Solo
            elif (nbTypes == 1) and ('solo' in types):
                bridgeType = 'solo'

            # c) Orphan
            elif (nbTypes == 1) and ('transduced' in types):
                bridgeType = 'orphan'

            # d) Unasigned type
            else:
                bridgeType = None
            
            ## Assess suppl. clusters bridge info consistency
            families = list(set([cluster.bridgeInfo['family'] for cluster in bridgeClusters['aligned'] if cluster.bridgeInfo['family']]))
            nbFamilies = len(families)

            srcIds = list(set([cluster.bridgeInfo['srcId'] for cluster in bridgeClusters['aligned'] if cluster.bridgeInfo['srcId']]))
            nbSrcIds = len(srcIds)

            ## a) Create bridge if consistent info:
            # 1. Bridge type determined
            # 2. Number of families not greater than 1
            # 3. Number of source elements not greater than 1
            if (bridgeType is not None) and (nbFamilies <= 1) and (nbSrcIds <= 1):
            
                ## Create dictionary containing bridge information
                bridgeInfo = {}
                bridgeInfo['bridgeType'] = bridgeType
                bridgeInfo['family'] = families[0]
                bridgeInfo['srcId'] = None if bridgeType == 'solo' else srcIds[0]
                bridgeInfo['bridgeLen'] = np.mean([cluster.bridgeInfo['bridgeLen'] for cluster in bridgeClusters['aligned']])
                bridgeInfo['bridgeSeq'] = None

                ## Create bridge
                self.bridge = BRIDGE(bridgeClusters['aligned'], bridgeInfo)

            ## b) Inconsistent info
            else:
                self.bridge = None

        ## 4. Apply read support filter to the bridge 
        ## Bridge available
        if self.bridge is not None:

            ## Compute percentage of metacluster supporting reads composing the bridge 
            percReads = float(self.bridge.nbReads()) / self.nbReadsTotal * 100

            ## Filter out bridge if:
            #    1) Bridge supported by < X reads OR
            #    2) Bridge supported by < X% of the total number of metacluster supporting reads 
            if (self.bridge.nbReads() < minReads) or (percReads < minPercReads):
                self.bridge = None


    def search4junctions(self, metaclustersBinDb, minReads, minPercReads):
        '''
        Use supplementary alignments to identify connections between the metacluster and any other BND metacluster. 

        Each of these connections will be a BND junction

        Input:
            1. metaclustersBinDb: bin database containing all the BND metaclusters identified in the sample
            2. minReads: minimum number of reads supporting the BND junction 
            3. minPercReads: minimum percentage of metacluster and partner metacluster reads supporting the BND junction

        Output:
            1. junctions: list of BND_junction objects
        '''
        ## 1. Search for candidate junction partners
        candidateJunctionPartners = []

        ## A) No bridge found
        if self.bridge is None:
            supplClusters = structures.dict2list(self.supplClusters)

            for cluster in supplClusters:

                # Candidate junction partner if adjacent
                if cluster.clusterIndex() == 0:
                    candidateJunctionPartners.append(cluster)

        ## Bridge found
        else:

            supportTypes = self.bridge.support_types()

            ## B) Bridge with unaligned suppl. cluster support
            if 'unaligned' in supportTypes:

                ## Second BND located at unaligned suppl. cluster
                cluster = self.bridge.return_unaligned_cluster()
                candidateJunctionPartners.append(cluster)

            ## C) Bridge with only aligned suppl. cluster support 
            else:                
                readIds = self.bridge.readIds()
                supplClusters = structures.dict2list(self.supplClusters)

                for cluster in supplClusters:

                    ## Assess if supplementary cluster is BND junction partner
                    partnerBool = cluster.is_junction_partner(readIds, minPercReads)

                    if partnerBool:
                        candidateJunctionPartners.append(cluster)

        ## 2. Select only those BND junction partners for which a metacluster is available
        candidateJunctions = []

        ## For each partner suppl. cluster
        for partner in candidateJunctionPartners:

            # Make sure there is a BND metacluster in the reference
            if partner.ref in metaclustersBinDb:

                ## Retrieve metacluster bin database on that reference
                binDbRef = metaclustersBinDb[partner.ref]

                ## Search for BND metaclusters at supplementary cluster breakpoint 
                buffer = 100
                partnerMetaclusters = binDbRef.collect_interval(partner.bkpPos() - buffer, partner.bkpPos() + buffer, ['BND'])

                ## Create BND junction candidate tuple for each metacluster - partner metacluster pair 
                for match in partnerMetaclusters:
                    metaPartner = match[0]
                    candidateJunctions.append((self, metaPartner, partner))

        ## 3. Apply support based filters and create BND junctions
        junctions = []

        ## For each junction
        for junction in candidateJunctions:
            meta, metaPartner, supplCluster = junction

            ## Compute % of metacluster and partner metacluster reads supporting BND junction 
            # Note: supplCluster is the suppl. cluster supporting the junction between 
            # the metaclusters
            nbReads = float(supplCluster.nbReads()[0])
            percReadsMeta = nbReads / meta.nbReadsTotal * 100
            percReadsPartner = nbReads / metaPartner.nbReadsTotal * 100 # Note: here % can be >100% (fix issue later by selecting suppl.cluster from the partner)

            ## Create BND junction if connection fulfills ALL these contitions: 
            #      1) Supported by >= minReads
            #      2) Percentage of reads supporting supplementary cluster >= minPercReads for both metacluster and partner
            if (nbReads >= minReads) and (percReadsMeta >= minPercReads) and (percReadsPartner >= minPercReads):
                                
                ## Create BND junction 
                junction = BND_junction(meta, metaPartner)
                
                ## Add bridge to the junction (None if no bridge found)
                junction.bridge = self.bridge

                ## Add junction to the list
                junctions.append(junction)

        return junctions
        

    def determine_INS_type(self, hits_genome, hits_splicing, hits_viral, repeatsDb, transducedDb, exonsDb, types2Search=['EXPANSION','DUP','INTERSPERSED','PSEUDOGENE']):
        '''
        Determine the type of insertion based on the alignments of the inserted sequence on the reference genome

        Input:
            1. hits_genome: PAF object containing inserted sequence alignments on the reference genome
            2. hits_splicing: list of bed entries containing inserted sequence alignments on the reference genome (splice-aware alignment). None if not hit found
            3. hits_viral: PAF object containing inserted sequence alignments on the viral database 
            4. repeatsDb: bin database containing annotated repeats in the reference. None if not available
            5. transducedDb: bin database containing regions transduced by source elements. None if not available
            6. exonsDb: bin database containing annotated exons. None if not available. 
            7. types2Search: Comma separated list containing insertion types to look for. Default: ['EXPANSION','DUP','INTERSPERSED','PSEUDOGENE']

        Output: Add INS type annotation to the attribute SV_features
        ''' 
        ## 0. Abort if consensus event not available 
        if self.consensusEvent is None and 'VIRUS' not in types2Search:
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0
            return 

        ## 1. Assess if input sequence corresponds to a repeat expansion
        if 'EXPANSION' in types2Search:
            is_EXPANSION, self.insertHits = self.is_expansion(hits_genome, repeatsDb)

            # Stop if insertion is a expansion
            if is_EXPANSION:
                return

        ## 2. Assess if input sequence corresponds to duplication 
        if 'DUP' in types2Search:
            is_DUP, self.insertHits = self.is_duplication(hits_genome, 100)

            # Stop if insertion is a duplication
            if is_DUP:
                return


        ## 3. Assess if input sequence corresponds to solo interspersed insertion or transduction
        # Note: return boolean as well specifying if interspersed or not
        if 'INTERSPERSED' in types2Search:
            is_INTERSPERSED, INS_features, self.insertHits = retrotransposons.is_interspersed_ins(self.consensusEvent.pick_insert(), hits_genome, repeatsDb, transducedDb)

            # Update metacluster with insertion features
            self.SV_features.update(INS_features) 

            # Stop if insertion is a interspersed insertion
            if is_INTERSPERSED:
                return    

        ## 4. Assess if input sequence corresponds to processed pseudogene insertion
        if 'PSEUDOGENE' in types2Search:
            is_PSEUDOGENE, outHits = self.is_processed_pseudogene(hits_splicing, exonsDb)

            # Stop if insertion is a processed pseudogene
            if is_PSEUDOGENE:
                return    

        ## 5. Assess if input sequence corresponds to a viral insertion
        if 'VIRUS' in types2Search:
            is_VIRUS, INS_features = self.is_VIRUS(hits_viral)

            # Stop if insertion is a viral insertion
            if is_VIRUS:
                # Update metacluster with insertion features
                self.SV_features.update(INS_features) 
                return is_VIRUS

    def is_VIRUS(self, PAF, side=None):
        '''
        Determine if input sequence corresponds to a viral insertion

        Input:
            1. PAF: PAF object containing input sequence alignments on the reference genome
            2. side: para los INS_noOverlap, saber que lado se esta analizando
        Output:
            1. VIRUS: Boolean specifying if inserted sequence corresponds to a virus (True) or not (False)
            2. INS_features: dictionary containing viral insertion features
        '''
        INS_features = {}
        identity = []
        specificIdentity = []
        

        print ('self.SV_type ' + str(self.SV_type))
        print ('self.beg ' + str(self.beg))

        ## 0a. Abort if no PAF available
        if not PAF:
            print ('Start not PAF')
            VIRUS = False
            INS_features['INS_TYPE'] = 'unknown'
            INS_features['PERC_RESOLVED'] = 0

            return VIRUS, INS_features

        ## 0b. Abort if no hit available
        if not PAF.alignments:
            print ('Start not PAF.alignments')
            VIRUS = False
            INS_features['INS_TYPE'] = 'unknown'
            INS_features['PERC_RESOLVED'] = 0

            return VIRUS, INS_features
        
        if self.SV_type == 'INS':
            # TODO: pedir un minimo
            #chain = PAF.chain(300, 20)
            chain = PAF.chain(300, 50)
            print ('PAFchain ' + str(chain.alignments))
            ## Identify it as a viral insertion if percentage of query covered is higher than 0.80
            percVirus = chain.perc_query_covered()
            # TODO: Poner in threshold!!
            if percVirus > 0:
                VIRUS = True
                INS_features['INS_TYPE'] = 'viral'
                INS_features['PERC_RESOLVED'] = percVirus
                for ali in PAF.alignments:
                    splitedtName = ali.tName.split('|')
                    identity.append(splitedtName[0])
                    if len(splitedtName) > 1:
                        specificIdentity.append(splitedtName[1])
                INS_features['IDENTITY'] = set(identity)
                if specificIdentity:
                    INS_features['SPECIDENTITY'] = specificIdentity
            else:
                VIRUS = False
                INS_features['INS_TYPE'] = 'unknown'
                INS_features['PERC_RESOLVED'] = 0
                INS_features['IDENTITY'] = None

        elif 'BND' in self.SV_type:
            VIRUS = True
            INS_features['INS_TYPE'] = 'viral'
            for ali in PAF.alignments:
                splitedtName = ali.tName.split('|')
                identity.append(splitedtName[0])
                if len(splitedtName) > 1:
                    specificIdentity.append(splitedtName[1])
            INS_features['IDENTITY'] = set(identity)
            if specificIdentity:
                INS_features['SPECIDENTITY'] = specificIdentity

        elif 'INS_noOverlap' in self.SV_type:
            VIRUS = True
            INS_features[side] = 'viral'
            INS_features['IDENTITY'] = 'viral'
            for ali in PAF.alignments:
                splitedtName = ali.tName.split('|')
                identity.append(splitedtName[0])
                if len(splitedtName) > 1:
                    specificIdentity.append(splitedtName[1])
            INS_features['IDENTITY_' + side] = set(identity)
            if specificIdentity:
                INS_features['SPECIDENTITY_' + side] = specificIdentity

        
        return VIRUS, INS_features

    def determine_INS_structure(self, consensusPath, transducedPath, transductionSearch, outDir):
        '''
        Infer inserted sequence structural features

        Input:
            1. consensusPath: path to fasta file containing retrotransposon consensus sequences
            2. transducedPath: path to fasta containing transduced sequences downstream of source elements
            3. transductionSearch: boolean specifying if transduction search is enabled (True) or not (False)
            4. outDir: output directory
    
        Output: 
            structure: dictionary containing insertion structural properties
        '''
        ##  Skip structure inference if consensus event not available
        if self.consensusEvent is None:
            return {}

        ## 1. Read fasta files 
        #  1.1 Consensus sequences
        consensus = formats.FASTA()
        consensus.read(consensusPath)

        #  1.2 Transduced regions
        if transductionSearch: 
            transduced = formats.FASTA()
            transduced.read(transducedPath)

        ## 2. Create fasta object containing database of sequences
        ## The database will contain the following sequences depending on the insertion type:
        ## - Solo      -> consensus sequences for the same family
        ## - Partnered -> consensus sequences for the same family
        #              -> corresponding transduced area
        ## - Orphan    -> corresponding transduced area

        ## Initialize fasta
        fasta = formats.FASTA()

        ## Add to the fasta subfamily consensus sequences for the corresponding family
        if self.SV_features['INS_TYPE'] in ['solo', 'partnered']:
            for seqId, seq in consensus.seqDict.items(): 
                family = seqId.split('|')[1]

                if family in self.SV_features['FAMILY']:
                    fasta.seqDict[seqId] = seq

        ## Add to the fasta transduced region or regions
        if self.SV_features['INS_TYPE'] in ['partnered', 'orphan']:
        
            for seqId, seq in transduced.seqDict.items(): 
                family, srcId = seqId.split('|')[1:3]

                if (family in self.SV_features['FAMILY']) and (srcId in self.SV_features['CYTOBAND']):
                    fasta.seqDict[seqId] = seq

        ##  Skip structure inference if empty database of sequences
        if not fasta.seqDict:
            return {}

        ## 3. Create fasta file
        fastaPath = outDir + '/reference_sequences.fa'
        fasta.write(fastaPath)

        ## 4. Index fasta file
        fileName = 'reference_sequences'  
        indexPath = alignment.index_minimap2(fastaPath, fileName, outDir)

        ## 5. Create fasta file containing consensus inserted sequence
        # Create fasta object
        FASTA = formats.FASTA()
        insert = self.consensusEvent.pick_insert()
        FASTA.seqDict['consensus_insert'] = insert

        # Write fasta
        insertPath = outDir + '/consensus_insert.fa'
        FASTA.write(insertPath)    

        ## 6. Structure inference        
        structure = retrotransposons.retrotransposon_structure(insertPath, indexPath, outDir)

        return structure

    def is_duplication(self, PAF, buffer):
        '''
        Determine if metacluster corresponds to a short tandem duplication
        
        Input:
            1. PAF: PAF object containing consensus inserted sequence alignments on the reference genome
            2. buffer: number of base pairs to extend metacluster interval when assessing overlap

        Output:
            1. DUP: boolean specifying if inserted sequence corresponds to a duplication (True) or not (False)
            2. HITS: PAF object containing inserted sequence alignments supporting a duplication. None if no hit found

        Update SV_features attribute with 'INS_TYPE' and 'PERC_RESOLVED' info

        Note: I need to modify the way PERC_RESOLVED. Now it´s not precise, do it based on the qBeg and qEnd
        alignment coordinates
        '''
        ## 0. Abort if no hit available
        if not PAF.alignments:
            DUP = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0

            return DUP, None   

        ## 1. Initialize
        totalPerc = 0
        HITS = formats.PAF()
        
        ## 2. Search hits matching metacluster interval
        # For each hit
        for hit in PAF.alignments: 

            ## Hit in the same ref as the metacluster
            if hit.tName == self.ref: 
                
                overlap, nbBp = gRanges.overlap(self.beg - buffer, self.end + buffer, hit.tBeg, hit.tEnd)
                perc = float(nbBp) / hit.qLen * 100

                ## Hit within metacluster interval
                if overlap:
                    totalPerc += perc
                    HITS.alignments.append(hit)

        ## set upper bound to 100
        if totalPerc > 100:
            totalPerc = 100 

        ## 3. Determine if duplication or not
        # a) Duplication
        if totalPerc >= 40:

            DUP = True
            self.SV_features['INS_TYPE'] = 'duplication'
            self.SV_features['PERC_RESOLVED'] = totalPerc

        # b) Not duplication
        else:

            DUP = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0

        return DUP, HITS

    def is_expansion(self, PAF, repeatsDb):
        '''
        Determine if metacluster corresponds to a repeat expansion
        
        Input:
            1. PAF: PAF object containing consensus inserted sequence alignments on the reference genome
            2. repeatsDb: bin database containing annotated repeats in the reference

        Output:
            1. EXPANSION: Boolean specifying if inserted sequence corresponds to an expansion (True) or not (False)
            2. HITS: PAF object containing inserted sequence alignments supporting an expansion. None if no hit found

        Update SV_features attribute with 'INS_TYPE', 'FAMILY', 'SUBFAMILY', 'PERC_RESOLVED' info

        Note: I need to modify the way PERC_RESOLVED. Now it´s not precise, do it based on the qBeg and qEnd
        alignment coordinates        
        '''

        ## 0. Abort if no hit or repeats database not available
        if (not PAF.alignments) or (repeatsDb is None):

            EXPANSION = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0

            return EXPANSION, None   
            
        ## 1. Initialize
        totalPerc = 0
        HITS = formats.PAF()

        ## 2. Assess if metacluster located over an annotated repeat sequence 
        # Make list of annotated repeat categories according to repeatmasker
        repeatTypes = ['Low_complexity', 'Simple_repeat', 'Satellite', 'telo', 'centr', 'acro']

        # Selecting those annotated repeats over the target region 
        repeatsFiltered = [ repeat for repeat in self.repeatAnnot if repeat['distance'] == 0 ]

        # Make list of annotated families at the target region
        targetFamilies = [ repeat['family'] for repeat in repeatsFiltered ]
        targetSubfamilies = [ repeat['subfamily'] for repeat in repeatsFiltered ]

        # A) Metacluster over annotated repeat
        if any(family in repeatTypes for family in targetFamilies):

            families = []
            subfamilies = []

            ### Search for inserted sequence hits on repeats of the same family
            # For each hit
            for hit in PAF.alignments: 

                ## Intersect hit coordinates with annotated repeats
                overlaps = annotation.annotate_interval(hit.tName, hit.tBeg, hit.tEnd, repeatsDb)
                
                ## For each repeat overlapping hit coordinates
                for overlap in overlaps:
                    event, nbBp, percBp, coord = overlap
                    
                    # Repeat belonging to the same family as the ones at metacluster interval
                    if event.optional['family'] in targetFamilies:
                        perc = float(nbBp) / hit.qLen * 100
                        totalPerc += perc
                        families.append(event.optional['family'])
                        subfamilies.append(event.optional['subfamily'])

                        if hit not in HITS.alignments:
                            HITS.alignments.append(hit)

            ## set upper bound to 100
            if totalPerc > 100:
                totalPerc = 100 

            ## Determine if expansion or not
            # a) Expansion
            if totalPerc >= 40:

                EXPANSION = True
                self.SV_features['INS_TYPE'] = 'expansion'
                self.SV_features['PERC_RESOLVED'] = totalPerc
                self.SV_features['FAMILY'] = list(set(families))
                self.SV_features['SUBFAMILY'] = list(set(subfamilies))

            # b) Not expansion
            else:

                EXPANSION = False
                self.SV_features['INS_TYPE'] = 'unknown'
                self.SV_features['PERC_RESOLVED'] = 0

        # B) Metacluster outside annotated repeat
        else:
            EXPANSION = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0

        return EXPANSION, HITS     

    def is_processed_pseudogene(self, hits, exonsDb):     
        '''
        Determine if metacluster corresponds to a processed pseudogene insertion
        
        Input:
            1. hits: list of bed entries corresponding to inserted sequence hits on the reference 
            2. exonsDb: bin database containing annotated exons in the reference. None if not available

        Output:
            1. PSEUDOGENE: Boolean specifying if inserted sequence corresponds to an expansion (True) or not (False)
            2. outHits: list of bed entries corresponding to hits matching annotated exons
        
        Update SV_features attribute with 'INS_TYPE', 'PERC_RESOLVED', 'NB_EXONS', 'SOURCE_GENE', 'POLYA', 'STRAND'

        Note: I need to modify the way PERC_RESOLVED. Now it´s not precise, do it based on the qBeg and qEnd
        alignment coordinates        
        '''

        ## 0. Abort if no hit available
        if (not hits) or (exonsDb is None):

            PSEUDOGENE = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0

            return PSEUDOGENE, None   

        ## 1. Search for polyA/T tails
        insert = self.consensusEvent.pick_insert()
        windowSize = 8
        maxWindowDist = 2
        minMonomerSize = 10
        minPurity = 80 

        # 1.1 Search for poly(A) monomers on the 3' end 
        targetSeq = insert[-50:]
        targetMonomer = 'A'
        monomers3end = sequences.find_monomers(targetSeq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

        ## 1.2 Search for polyT on the 5' end 
        targetSeq = insert[:50]
        targetMonomer = 'T'
        monomers5end = sequences.find_monomers(targetSeq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

        if monomers3end and monomers5end:
            polyA = False
            strand = None

        elif monomers3end:
            polyA = True
            strand = '+'

        elif monomers5end:
            polyA = True
            strand = '-'

        else:
            polyA = False
            strand = None

        ## 2. Intersect each hit with the annotated exons
        outHits = []
        nbResolvedBp = 0
        nbExons = 0

        # For each hit
        for hit in hits: 

            hit.annot = {}

            ## Do intersection
            overlaps = annotation.annotate_interval(hit.ref, hit.beg, hit.end, exonsDb)

            ## Hit intersect with exon
            if overlaps:
                longestOverlap = overlaps[0] # Select exon with longest overlap
                hit.annot['EXON'] = longestOverlap[0] 
                
                outHits.append(hit)
                nbResolvedBp += longestOverlap[1]
                nbExons += 1

        ## 3. Compute percentage of inserted sequence matching on annotated exons
        percResolved = float(nbResolvedBp) / self.consensusEvent.length * 100

        # set upper bound to 100
        if percResolved > 100:
            percResolved = 100 

        ## 4. Determine if pseudogene or not
        # a) Pseudogene (enough % of sequence resolved + polyA)
        if (percResolved >= 60) and (polyA):

            PSEUDOGENE = True
            self.SV_features['INS_TYPE'] = 'pseudogene'
            self.SV_features['PERC_RESOLVED'] = percResolved
            self.SV_features['NB_EXONS'] = nbExons
            self.SV_features['SOURCE_GENE'] = list(set([hit.annot['EXON'].optional['geneName'] for hit in outHits]))
            self.SV_features['POLYA'] = polyA
            self.SV_features['STRAND'] = strand

        # b) Not pseudogene
        else:
            PSEUDOGENE = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = percResolved        

        return PSEUDOGENE, outHits

    def setElement(self):
        self.element = None
        for event in self.events:
            if hasattr(event, 'element'):
                self.element = event.element
                break
        return self.element
    
class BRIDGE():
    '''
    Rearrangement bridge class
    '''
    def __init__(self, supplClusters, info):

        self.supplClusters = supplClusters
        self.bridgeType = info['bridgeType']
        self.family = info['family']
        self.srcId = info['srcId']
        self.bridgeSeq = info['bridgeSeq']
        self.bridgeLen = info['bridgeLen']

    def add(self, supplClusters):
        '''
        Incorporate additional supplementary clusters into the bridge

        Input:
            1. supplClusters: list of supplementary clusters
        '''
        self.supplClusters.extend(supplClusters)
        
    def nbClusters(self):
        '''
        Return the number of supplementary clusters composing the bridge
        '''
        return len(self.supplClusters)

    def readIds(self):
        '''
        Return a list containing bridge supporting reads 
        '''
        ## Collect read ids supporting the suppl. clusters composing the bridge
        readIdsNested = [cluster.nbReads()[1] for cluster in self.supplClusters]
        readIds = list(itertools.chain.from_iterable(readIdsNested))

        ## Remove redundant read ids
        readIds = list(set(readIds))

        return readIds

    def nbReads(self):
        '''
        Return the number of reads supporting a rearrangement bridge
        '''

        ## Collect read ids supporting the suppl. clusters composing the bridge
        readIdsNested = [cluster.nbReads()[1] for cluster in self.supplClusters]
        readIds = list(itertools.chain.from_iterable(readIdsNested))

        ## Remove redundant read ids
        readIds = list(set(readIds))

        ## Compute number of distinct reads supporting the bridge
        nbReads = len(readIds)

        return nbReads
            
    def support_types(self):
        '''
        Return list of supplementary clusters types supporting the bridge (aligned or unaligned)
        '''
        supportTypes = []

        for cluster in self.supplClusters:

            if cluster.bridgeInfo['supportType'] not in supportTypes:
                supportTypes.append(cluster.bridgeInfo['supportType'])

        return supportTypes

    def return_unaligned_cluster(self):
        '''
        Return supplementary cluster supporting unaligned bridge
        '''
        targetCluster = None

        for cluster in self.supplClusters:
        
            if cluster.bridgeInfo['supportType'] == 'unaligned':
                targetCluster = cluster
        
        return targetCluster


class BND_junction():
    '''
    Rearrangement bridge class
    '''
    def __init__(self, metaclusterA, metaclusterB):

        self.metaclusterA = metaclusterA
        self.metaclusterB = metaclusterB

        ## Initialize bridge information
        self.bridge = None

        ## Consensus information
        self.junctionConsSeq = None
        self.bkpsConsSupport = None # Boolean telling if metaclsuters bkps are supported by consensus realignment with buffr of 100 bp

        self.identity = None

    def junctionCoord(self):
        '''
        Return a string containing junction coordinates with the following format: 
            refA:bkpCoordA-refB:bkpCoordB
        '''
        coord = self.metaclusterA.ref + ':' + str(self.metaclusterA.bkpPos) + '-' + self.metaclusterB.ref + ':' + str(self.metaclusterB.bkpPos)
        return coord

    def junctionCoord_rev(self):
        '''
        Return a string containing junction coordinates with the following format: 
            refB:bkpCoordB-refA:bkpCoordA        
        '''
        coord = self.metaclusterB.ref + ':' + str(self.metaclusterB.bkpPos) + '-' + self.metaclusterA.ref + ':' + str(self.metaclusterA.bkpPos)
        return coord

    def junctionType(self):
        '''
        Determine BND junctions type.

        Output:
            1. junctionType: intrachromosomal or interchromosomal (for now... later make code able to differenciate tanden dup, inversions, ...)
        '''

        # 1) Intrachromosomal
        if self.metaclusterA.ref == self.metaclusterB.ref:
            junctionType = 'intrachromosomal'

        # 2) Interchromosomal 
        else:
            junctionType = 'interchromosomal'

        return junctionType             
        
    def supportingReads(self):
        '''
        Compute the number of reads supporting the BND junction

        Output:
            1. nbTotal: Total number of BND junction supporting reads
            2. nbTumour: Number of BND junction supporting reads in the tumour 
            3. nbNormal: Number of BND junction supporting reads in the normal
        '''
        ## Total number of reads
        nbTotal = len(set(self.metaclusterA.reads + self.metaclusterB.reads))

        ## Number of reads in the tumour

        # NOTE 2020: change 2020 check it!!!
        if self.metaclusterA.readsTumour != None:


            nbTumour = len(set(self.metaclusterA.readsTumour + self.metaclusterB.readsTumour))

            ## Number of reads in the normal
            nbNormal = len(set(self.metaclusterA.readsNormal + self.metaclusterB.readsNormal))
        
        else:
            nbTumour = None
            nbNormal = None


        return nbTotal, nbTumour, nbNormal
    
    def extractSupportingRead(self):
        '''
        Method for extracting the events that completely span the BND_junction 
        Output: 
            - junctionsList:
                - field 0: list of metaclusters that contain the same read
                - field 1: times that the read appears in a different events
                - field 2: number of supplementary alignments of the read.
        '''

        readNamesJunction={}
        '''
        - Dictionary:
        - keys: read names supporting the BND_junction
        - values: list:
            - field 0: list containing event objects corresponding to read name (only metaclsuterA and metaclsuterB events are kept, no bridge events)
            - field 1: number of times that this read name is repeated in the BND_junction components
            - field 2: number of read supplementary alignments 
        '''

        # 1. Make a list of readNames of events of metaclusterA
        for event in self.metaclusterA.events:
            try:
                readNamesJunction[event.readName]=[[event], 1, event.supplAlignment.count(';')]
            except AttributeError:
                readNamesJunction[event.readName]=[[event], 1, 0]

        # 2. Add to the list of readNames those from events of metaclusterB
        for event in self.metaclusterB.events:
            try:
                readNamesJunction[event.readName][0].append(event)
                readNamesJunction[event.readName][1] = readNamesJunction[event.readName][1] + 1
            except KeyError:
                try:
                    readNamesJunction[event.readName]=[[event], 1, event.supplAlignment.count(';')]
                except AttributeError:
                    readNamesJunction[event.readName]=[[event], 1, 0]

        if self.bridge:
            # 3. Add to the list of readNames those from events of bridge
            for SUPPLEMENTARY_cluster in self.bridge.supplClusters:
                for event in SUPPLEMENTARY_cluster.events:
                    try:
                        readNamesJunction[event.readName][1] = readNamesJunction[event.readName][1] + 1
                    except KeyError:
                        readNamesJunction[event.readName]=[event, 1, event.supplAlignment.count(';')]
        
        # TODO: delete readNamesJunction

        # 4. Make list from dictionary, excluding readNames.
        junctionsList = []
        for junctionsLists in readNamesJunction.values():
            junctionsList.append(junctionsLists)
        junctionsList = sorted(junctionsList, key = lambda x: (-x[1], x[2]))

        return junctionsList

    def get_consensus_BNDInterval(self, polishedFastaEntireSequence, junctionInterval, outDir):
        '''
        From a BND_junction consensus fasta file, get only the sequence containing bkps +- 1000
        Input:
            - polishedFastaEntireSequence: consensus fasta file of the entire sequence
            - junctionInterval: bkps of BND_junction
        Output:
            - polishedFastaInterval: Path to the consensus BND_junction fasta file
            - polishedFastaIntervalObj: object of the consensus BND_junction fasta file
        '''

        # Pick from polished fasta only the region involving the junction (using a buffer of 1000bp)
        polishedFastaIntervalDict = {}
        polishedFastaObj = formats.FASTA()
        polishedFastaObj.read(polishedFastaEntireSequence)
        for key, value in polishedFastaObj.seqDict.items():
            beg = min(junctionInterval)-5000 if min(junctionInterval) > 5000 else 0
            polishedFastaIntervalDict[key]=value[beg:max(junctionInterval)+5000]

        polishedFastaInterval = polishedFastaEntireSequence + 'fasta' 
        polishedFastaIntervalObj = formats.FASTA()
        polishedFastaIntervalObj.seqDict = polishedFastaIntervalDict
        polishedFastaIntervalObj.write(polishedFastaInterval)

        return polishedFastaInterval, polishedFastaIntervalObj
    
    def get_BNDjunction_chain(self, polishedFastaInterval, target, viralDb, outDir):
        '''
        Perform local aligments of BND_junction and get a PAF chain
        Input: 
            - polishedFastaInterval: Path to the consensus BND_junction fasta file
            - target: Path to targeted reference regions
            - outDir
        Output:
            - junctionPAFChain: consensus BND_junction PAF chain (None no chain made)
        '''
        junctionPAFChain = None
        if self.bridge:
            # Merge target and tr fasta in one:
            # NOTE 2020: SILENCE 2020 DESILENCE
            # TODO: Make the difference bewteen partnered, orphan and repeat
            '''
            if self.bridge.bridgeType == 'partnered' or self.bridge.bridgeType == 'orphan' or self.bridge.bridgeType == 'repeat':
                targetFasta = formats.FASTA()
                targetFasta.read(target)

                refSeqsIndex = outDir + '/reference_sequences.fa'
                refSeqsIndexFasta = formats.FASTA()
                refSeqsIndexFasta.read(refSeqsIndex)

                trueFasta = formats.merge_FASTA([targetFasta, refSeqsIndexFasta])
                trueFastaPath = outDir + 'trueFasta.fa'
                trueFasta.write(trueFastaPath)

                junctionPAFPath = alignment.alignment_minimap2(polishedFastaInterval, trueFastaPath, 'alljunctions', 1, outDir)
                if os.path.getsize(junctionPAFPath) > 0:

                    # Make junction chain
                    junctionPAF = formats.PAF()
                    junctionPAF.read(junctionPAFPath)

                    # TODO: ensure that these parameters are ok!!!
                    junctionPAFChain = junctionPAF.chain(50,30)
            '''
            # NOTE 2020: Change 2020
            if 'viral' in self.bridge.bridgeType:
                targetFasta = formats.FASTA()
                targetFasta.read(target)

                # TODO: change db path
                refSeqsIndex = viralDb
                refSeqsIndexFasta = formats.FASTA()
                refSeqsIndexFasta.read(refSeqsIndex)

                trueFasta = formats.merge_FASTA([targetFasta, refSeqsIndexFasta])
                trueFastaPath = outDir + 'trueFasta.fa'
                trueFasta.write(trueFastaPath)

                junctionPAFPath = alignment.alignment_minimap2(polishedFastaInterval, trueFastaPath, 'alljunctions', 1, outDir)
                if os.path.getsize(junctionPAFPath) > 0:

                    # Make junction chain
                    junctionPAF = formats.PAF()
                    junctionPAF.read(junctionPAFPath)

                    # TODO: ensure that these parameters are ok!!!
                    junctionPAFChain = junctionPAF.chain(100,30)

        else:
            junctionPAFPath = alignment.alignment_minimap2(polishedFastaInterval, target, 'alljunctions', 1, outDir)
            if os.path.getsize(junctionPAFPath) > 0:
                # Make junction chain
                junctionPAF = formats.PAF()
                junctionPAF.read(junctionPAFPath)

                # TODO: ensure that these parameters are ok!!!
                junctionPAFChain = junctionPAF.chain(50,30)

        return junctionPAFChain

    def get_consensus_BNDSeq(self, junctionPAFChain, polishedFastaIntervalObj):
        '''
        Get junction consensus sequence(junctionConsSeq) and names and coordinates of bridge sequences
        Input:
            - junctionPAFChain
            - polishedFastaIntervalObj: object of the consensus BND_junction fasta file
        Output:
            - aligTNames: list of names of bridge sequences
            - aligCoordinates: list of read coordinates of bridge sequences
        '''
        consensusSequence = ""
        aligTNames = []
        aligCoordinates = []
        for alig2 in junctionPAFChain.alignments:
            # Mark consensus sequence and save it
            # Read polished fasta:
            polishedSequence = list(polishedFastaIntervalObj.seqDict.values())[0]

            polishedSequenceSeg = polishedSequence[min([alig2.qBeg,alig2.qEnd]):max([alig2.qBeg,alig2.qEnd])]
            if consensusSequence == "":
                consensusSequence = polishedSequenceSeg
            else:
                consensusSequence = consensusSequence + '|' + polishedSequenceSeg
            # TODO: This is a temp fix for unclass that has an space on its name, but it should be fixed in another way
            if '|' in alig2.tName or 'unclass' in alig2.tName or 'environ' in alig2.tName:

                aligTNames.append(alig2.tName)
                aligCoordinates.append(alig2.qBeg)
                aligCoordinates.append(alig2.qEnd)
        
        # Save consSeq +-1000 indicating where is the bkp ("|") and chain as attributes de BND_junction
        self.junctionConsSeq = consensusSequence

        return aligTNames, aligCoordinates

    def analyse_BNDjunction_bridge(self, aligTNames, aligCoordinates, polishedFastaIntervalObj, viralDb, outDir):
        '''
        Realing BND_junction bridge and analyse its sequence.
        - Input:
            - aligTNames: list of names of bridge sequences
            - aligCoordinates: list of read coordinates of bridge sequences
            - polishedFastaIntervalObj: object of the consensus BND_junction fasta file
            - outDir
        - Output:
            This function only modifies bridge features:
        '''

        # TODO: Check if we are loosing something due to this second condition:
        if self.bridge and len (aligTNames) > 0:
            # Pick from polishedFastaInterval only the region involving the bridge
            polishedFastaBridgeDict = {}
            polishedFastaBridgeObj = formats.FASTA()
            for key, value in polishedFastaIntervalObj.seqDict.items():
                polishedFastaBridgeDict[key]=value[min(aligCoordinates):max(aligCoordinates)]
            polishedFastaBridgePath = outDir + 'polishedFastaBridge.fa'
            polishedFastaBridgeObj.seqDict = polishedFastaBridgeDict
            polishedFastaBridgeObj.write(polishedFastaBridgePath)

            if any('consensus' in x for x in aligTNames) or any('transduced' in x for x in aligTNames):
                # Look for retrotransposon_structure
                # NOTE (EVA): This function was run before (in supports_unaligned_bridge) and I think it's quite redundant
                refSeqsIndex = outDir + '/reference_sequences.fa'
                structure = retrotransposons.retrotransposon_structure(polishedFastaBridgePath, refSeqsIndex, outDir)


                # a) Resolved structure
                if ('INS_TYPE' in structure) and (structure['INS_TYPE'] is not 'unknown') and ('PERC_RESOLVED' in structure) and (structure['PERC_RESOLVED'] >= 60):
                    self.bridge.bridgeType = structure['INS_TYPE'] 
                    # TODO: put consensus sequence
                    #self.bridge.bridgeInfo['bridgeSeq'] = bridgeSeq
                    #self.bridge.bridgeInfo['bridgeLen'] = bridgeLen
                    self.bridge.family = ','.join(structure['FAMILY']) 
                    self.bridge.srcId = ','.join(structure['CYTOBAND']) if ('CYTOBAND' in structure and structure['CYTOBAND']) else None
        
            # For the moment, only Rt and viruses have this structure:
            else:
                # Look for viral structure
                # TODO: change db path
                # NOTE (EVA): This function was run before (in supports_unaligned_bridge) and I think it's quite redundant
                structure = virus.virus_structure(polishedFastaBridgePath, viralDb, outDir)

                if ('INS_TYPE' in structure) and (structure['INS_TYPE'] is not 'unknown') and ('PERC_RESOLVED' in structure) and (structure['PERC_RESOLVED'] >= 60):
                    # b) Unresolved structure
                    ## TODO: VIRUSES STRUCTURE!!!
                    self.bridge.bridgeType = structure['INS_TYPE'] 
                    # TODO: put consensus sequence
                    #self.bridge.bridgeInfo['bridgeSeq'] = bridgeSeq
                    #self.bridge.bridgeInfo['bridgeLen'] = bridgeLen
                    self.bridge.family = ','.join(structure['FAMILY']) 
                    self.bridge.srcId = ','.join(structure['CYTOBAND']) if ('CYTOBAND' in structure and structure['CYTOBAND']) else None
                else:
                    self.bridge.bridgeType = 'unknown'
                    self.bridge.family = None
                    self.bridge.srcId = None

    def consensusBNDjunction_support_Bkps(self, junctionPAFChain):
        '''
        Check if consensus BND_junction realignment against reference supports previuos bkps.
        Input:
            - junctionPAFChain
        Output:
            - Boolean: True if bkps are supported, False otherwise
        '''

        # NOTE 2020: New 2020

        self.bkpsConsSupport = False

        if junctionPAFChain.alignments:

            # NOTE 2020: New 2020
            # TODO 2020: Put in a good way
            if ':' in junctionPAFChain.alignments[0].tName and ':' in junctionPAFChain.alignments[-1].tName:

                print ('junctionPAFChain.alignments[0].tName ' + str(junctionPAFChain.alignments[0].tName))
                print ('junctionPAFChain.alignments[-1].tName ' + str(junctionPAFChain.alignments[-1].tName))
                # Check if consensus alignment supports junction structure.
                consensusFirsttName = junctionPAFChain.alignments[0].tName.split(':')[1].split('-')[0]
                consensusLasttName = junctionPAFChain.alignments[-1].tName.split(':')[1].split('-')[0]

                consensusFirstStrand = junctionPAFChain.alignments[0].strand
                consensusLastStrand = junctionPAFChain.alignments[-1].strand

                if consensusFirstStrand == '+':
                    consensusBeg = int(consensusFirsttName) + junctionPAFChain.alignments[0].tEnd
                else:
                    consensusBeg = int(consensusFirsttName) + junctionPAFChain.alignments[0].tBeg

                if consensusLastStrand == '+':
                    consensusEnd = int(consensusLasttName) + junctionPAFChain.alignments[-1].tBeg
                else:
                    consensusEnd = int(consensusLasttName) + junctionPAFChain.alignments[-1].tEnd

                # Sort in order to compare always the same:
                metaclustersBkps = [self.metaclusterA.bkpPos, self.metaclusterB.bkpPos]
                consensusBkps = [consensusBeg,consensusEnd]
                metaclustersBkps.sort()
                consensusBkps.sort()

                # Two ifs: Check if always consensusBeg is always bkpA and viceversa
                if abs(metaclustersBkps[0] - consensusBkps[0]) > 1500 or abs(metaclustersBkps[1] - consensusBkps[1]) > 1500:
                    print ('Consensus aligment dont support junction structure')
                    self.bkpsConsSupport = False
                else:
                    print ('Consensus aligment support junction structure')
                    self.bkpsConsSupport = True

def metacluster_mate_suppl(discordants, leftClippings, rightClippings, minReads, refLengths):
    '''
    Group discordant and clipping clusters into metaclusters

    Input:
        1. discordants: list of discordant cluster objects
        2. leftClippings: list of left clipping cluster objects
        3. rightClippings: list of right clipping cluster objects
        4. minReads: minimum number of reads
        5. refLengths: dictionary containing references as keys and their lengths as values

    Output:
        1. filteredMeta: list of metaclusters
    '''
    ## 1. Create list of discordant mate clusters
    mateClusters = [discordant.create_matesCluster() for discordant in discordants]
    
    ## 2. Create list of supplementary clusters
    supplClustersLeft = [clipping.supplCluster for clipping in leftClippings]
    supplClustersRight = [clipping.supplCluster for clipping in rightClippings]
    supplClusters = supplClustersLeft + supplClustersRight
    
    ## 2.1. Fill supplementary cluster orientation attribute
    [supplCluster.clippOrientation() for supplCluster in supplClusters]
    
    ## 3. Organize discordant mate and suppl. clusters into a dictionary
    mateDict = events.events2nestedDict(mateClusters, 'DISCORDANT')
    supplDict = events.events2nestedDict(supplClusters, 'SUPPLEMENTARY')
    clustersDict = events.mergeNestedDict(mateDict, supplDict)

    ## 4. Organize discordant mate and suppl. clusters into a bin database
    wgBinDb = structures.create_bin_database(refLengths, clustersDict)

    ## 5. Perform metaclustering
    allMeta = []

    for binDb in wgBinDb.values():
        
        meta = clustering.reciprocal_overlap_clustering(binDb, 1, 1, binDb.eventTypes, 500, 'META')
        allMeta = allMeta + meta

    ## 6. Filter metaclusters based on read support
    filteredMeta = []

    for meta in allMeta:

        if meta.supportingReads()[0] >= minReads:
            filteredMeta.append(meta)
    
    return filteredMeta


def discCluster_from_readNames(readNames, ref, beg, end, bam, sample):
    '''
    Create discordant cluster from indicated reads in the region ref:beg-end. Useful to create cluster of mates having all the info.
    
    Input:
    1. readNames: List of read names
    2. ref
    3. beg
    4. end
    5. bam: Bam file
    6. sample: Sample type ['TUMOUR', 'NORMAL', None]
    
    Output: 
    1. newCluster: Cluster of discordants 
    '''
    DISCORDANTS = []
    
    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, beg, end)
    
    # For each read alignment
    for alignmentObj in iterator:
        
        # if end coordinate is None, continue. 
        # returns None if read is unmapped or no cigar alignment present
        if alignmentObj.reference_end == None:
            continue
        
        if alignmentObj.query_name in readNames:
            
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

            DISCORDANT = events.DISCORDANT(alignmentObj.reference_name, alignmentObj.reference_start, alignmentObj.reference_end, orientation, pair, alignmentObj.query_name, alignmentObj, sample, alignmentObj.is_duplicate)
            DISCORDANTS.append(DISCORDANT)
        
    if DISCORDANTS:        
        newCluster = create_cluster(DISCORDANTS, 'DISCORDANT')
    else:
        newCluster = None
    
    return(newCluster)