'''
Module 'filters' - Contains functions for filtering clusters
'''
## External
import pysam
import statistics

## Internal
from GAPI import gRanges
from GAPI import bamtools
from GAPI import stats
from GAPI import structures
from GAPI import clusters
from GAPI import events


###############
## FUNCTIONS ##
###############

def filter_clusters(clusters, filters2Apply, bam, normalBam, confDict, clusterType):
    '''
    Function to apply filters all discordant clusters. 

    Input:
        1. clusters: list of clusters
        2. filters2Apply: list containing the filters to apply 
        3. bam: path to bam file. None if not needed
        4. normalBam: path to matched normal bam file. None if not needed        
        5. confDict
        6. clusterType: 'CLIPPING', 'DISCORDANT', 'META'

    Output:
        1. clustersPass: List of clusters passing all the filters
    '''
    clustersPass = []
    
    # For each cluster
    for cluster in clusters:

        ## Apply filters
        if clusterType == 'DISCORDANT':
            failedFilters = filter_discordant(cluster, filters2Apply, bam, normalBam, confDict)
            
        elif clusterType == 'CLIPPING':
            failedFilters = filter_clipping(cluster, filters2Apply, confDict)
        
        elif clusterType == 'META':
            failedFilters = filter_metacluster(cluster, filters2Apply, confDict, bam, normalBam, mode='SR')
                  
        # Cluster pass all the filters
        if not failedFilters: 
            clustersPass.append(cluster)
            
        else:
            print("Discarded", clusterType, "cluster")
            print(cluster.ref, cluster.beg, cluster.end)
            print(failedFilters)
            print([event.readName for event in cluster.events])
    
    return clustersPass

def filter_clipping(clipping, filters2Apply, confDict):
    '''
    Apply selected filters to a clipping cluster provided as input

    Input:
        1. clipping: clipping cluster object
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        5. confDict:

    Output:
        1. failedFilters -> list containing those filters that the clipping cluster didn´t pass.
    '''        
    failedFilters = []

    ## 1. FILTER 1: Minimum number of reads per cluster
    if 'MIN-NBREADS' in filters2Apply: 
        if not filter_min_nb_reads(clipping, confDict['minNbCLIPPING'], confDict['minNormalReads']):
            failedFilters.append('MIN-NBREADS')
        
    ## 2. FILTER 2: Filter out those clusters with suppl outside target reference ##
    if 'SUPPL-REF' in filters2Apply: 

        if not filter_clipping_suppl_ref(clipping, confDict['targetRefs']):
            failedFilters.append('SUPPL-REF')

    ## 3. FILTER 3: Filter out those clusters whose supplementary alignments map over any source element downstream region 
    if 'SUPPL-SRC' in filters2Apply:

        if not filter_clipping_suppl_position(clipping, confDict['rangesDict'], 10000):
            failedFilters.append('SUPPL-SRC')
    
    ## 4. FILTER 4: Average mapping quality of supplementary alignments 
    if 'SUPPL-MAPQ' in filters2Apply:
    
        if not filter_suppl_MAPQ(clipping, 10):
            failedFilters.append('SUPPL-MAPQ')

    ## 5. FILTER 5: filter out clusters formed by tumour and normal reads. Discard germline variation
    if 'GERMLINE' in filters2Apply:

        if not filter_germline(clipping, confDict['minNormalReads']):
            failedFilters.append('GERMLINE')

    ## 6. FILTER 6: Filter out clusters based on duplicate percentage (Ex: 40%) 
    if 'READ-DUP' in filters2Apply:

        if not filter_highDup_clusters(clipping, 75):
            failedFilters.append('READ-DUP')
    
    ## 7. FILTER 7: Filter out clusters based on cluster coordinates ##
    if 'CLUSTER-RANGE' in filters2Apply:
        
        if not filter_clusterRange_clipping(clipping):
            failedFilters.append('CLUSTER-RANGE')

    return failedFilters

def filter_discordant(discordant, filters2Apply, bam, normalBam, confDict):
    '''
    Apply selected filters to a discordant cluster provided as input

    Input:
        1. discordant: discordant read pair cluster object
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. bam: path to bam file. None if not needed
        4. normalBam: path to matched normal bam file. None if not needed
        5. confDict

    Output:
        1. failedFilters -> list containing those filters that the discordant cluster doesn't pass.
    '''        
    failedFilters = []

    ## 1. FILTER 1: Minimum number of reads per cluster
    if 'MIN-NBREADS' in filters2Apply: 
        if not filter_min_nb_reads(discordant, confDict['minNbDISCORDANT'], confDict['minNormalReads']):
            failedFilters.append('MIN-NBREADS')

    ## 2. FILTER 2: Filter out those clusters with mates not over NOT target reference ##
    if 'MATE-REF' in filters2Apply: 

        if not filter_discordant_mate_ref(discordant, confDict['targetRefs']):
            failedFilters.append('MATE-REF')

    ## 3. FILTER 3: Filter out those clusters whose mates aligns over any source element downstream region ##
    if 'MATE-SRC' in filters2Apply:

        if not filter_discordant_mate_position(discordant, confDict['rangesDict'], 10000):
            failedFilters.append('MATE-SRC')
        
    ## 4. FILTER 4: Filter out clusters based on average MAPQ for mate alignments ##
    if 'MATE-MAPQ' in filters2Apply:
    
        if not filter_discordant_mate_MAPQ(discordant, 10, bam, normalBam):
            failedFilters.append('MATE-MAPQ')
        
    ## 5. FILTER 5: filter out clusters formed by tumour and normal reads. Discard germline variation (TILL HERE) 
    if 'GERMLINE' in filters2Apply:

        if not filter_germline(discordant, confDict['minNormalReads']):
            failedFilters.append('GERMLINE')
            
    ## 6. FILTER 6: Filter out clusters in unspecific regions ##
    if 'UNSPECIFIC' in filters2Apply:

        if not filter_discordant_mate_unspecific(discordant, 0.95, bam):
            failedFilters.append('UNSPECIFIC')

    ## 7. FILTER 7: Filter out clusters based on duplicate percentage (Ex: 50%) ##
    if 'READ-DUP' in filters2Apply:

        if not filter_highDup_clusters(discordant, 75):
            failedFilters.append('READ-DUP')
            
    ## 8. FILTER 8: Filter out clusters based on mates beg coordinates ##
    if 'CLUSTER-RANGE' in filters2Apply:
        
        if not filter_clusterRange_discordant(discordant):
            failedFilters.append('CLUSTER-RANGE')
    
    return failedFilters


def filter_metaclusters_sonia(metaclusters, filters2Apply, confDict, bam, normalBam):
    '''
    Function to apply filters all metaclusters. 

    Input:
        1. metaclustersDict: dictionary with the following structure: keys -> SV_type, value -> list of metaclusters corresponding to this SV_type.
        2. filters2Apply: dictionary containing lists as values list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. confDict
        4. bam

    Output:
        1. metaclustersPassDict: Dictionary with same structure as the input one, containing those metaclusters that passed all the filters.
        2. metaclustersFailDict: Dictionary with same structure as the input one, containig those metaclusters that failed one or more filters.
    '''

    ## For each metacluster
    for metacluster in metaclusters:
        
        # Apply filters
        metacluster.failedFilters = filter_metacluster(metacluster, filters2Apply, confDict, bam, normalBam)


def filter_metaclusters(metaclustersDict, filters2Apply, confDict, mode='SR'):
    '''
    Function to apply filters to a set of metaclusters organized in a dictionary

    Input:
        1. metaclustersDict: dictionary with the following structure: keys -> SV_type, value -> list of metaclusters corresponding to this SV_type.
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. confDict

    Output:
        1. metaclustersPassDict: Dictionary with same structure as the input one, containing those metaclusters that passed all the filters.
        2. metaclustersFailDict: Dictionary with same structure as the input one, containig those metaclusters that failed one or more filters.
    '''

    metaclustersPassDict = {}
    metaclustersFailDict = {}

    ## For each type of SV
    for SV_type, metaclusters in metaclustersDict.items():

        ## 1. Make list with the indexes of the metaclusters do not passing some filter
        filteredIndexes = []

        ## For each metacluster
        for index, metacluster in enumerate(metaclusters):

            ## Apply filters
            metacluster.failedFilters = filter_metacluster(metacluster, filters2Apply, confDict, None, mode=mode)

            # Metacluster fails some filter
            if metacluster.failedFilters:
                filteredIndexes.append(index)

        ## 2. Divide metaclusters in those passing and failing filtering
        for index, metacluster in enumerate(metaclusters):
            
            ## a) Failing some filter
            if index in filteredIndexes:

                ## Initialize list
                if SV_type not in metaclustersFailDict:
                    metaclustersFailDict[SV_type] = []
                
                metaclustersFailDict[SV_type].append(metacluster)

            ## b) Passing all the filters
            else:

                ## Initialize list
                if SV_type not in metaclustersPassDict:
                    metaclustersPassDict[SV_type] = []

                metaclustersPassDict[SV_type].append(metacluster)

    return metaclustersPassDict, metaclustersFailDict

def filter_metacluster(metacluster, filters2Apply, confDict, bam, normalBam = None, mode='SR'):
    '''
    Apply selected filters to one metacluster.

    Input:
        1. metacluster: metacluster object
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. confDict
        4. bam  
        5. normalBam

    Output:
        1. failedFilters -> list containing those filters that the metacluster doesn't pass.
    '''        
    failedFilters = []

    ## 1. FILTER 1: Minimum number of reads per cluster
    if 'MIN-NBREADS' in filters2Apply: 
        if not filter_min_nb_reads(metacluster, confDict['minReads'], confDict['minNormalReads']):
            failedFilters.append('MIN-NBREADS')

    ## 2. FILTER 2: Maximum number of reads per cluster
    if 'MAX-NBREADS' in filters2Apply: 
        if not filter_max_nb_reads(metacluster, confDict['maxClusterSize']):
            failedFilters.append('MAX-NBREADS')

    ## 3. FILTER 3: Maximum Coefficient of Variance per cluster
    if ('CV' in filters2Apply) and ('INS' in metacluster.subclusters): 
        if not filter_max_cv(metacluster, confDict['maxClusterCV']):
            failedFilters.append('CV')

    ## 4. FILTER 4: Whether a metacluster has a SV_type assigned or not
    if 'SV-TYPE' in filters2Apply: 
        if not filter_SV_type(metacluster, confDict['targetSV']):
            failedFilters.append('SV-TYPE')

    ## 5. FILTER 5: Minimum percentage of inserted sequence resolved
    if ('PERC-RESOLVED' in filters2Apply) and (metacluster.SV_type == 'INS') and ('PERC_RESOLVED' in metacluster.SV_features): 

        if not filter_perc_resolved(metacluster, confDict['minPercResolved']):
            failedFilters.append('PERC-RESOLVED')
            
    ## 6. FILTER 6: Area mapping quality
    ## 7. FILTER 7: Area clipping SMS
    if 'AREAMAPQ' in filters2Apply or 'AREASMS' in filters2Apply:
        areaMeta = area(metacluster, confDict, bam)
        if 'AREAMAPQ' in filters2Apply and not areaMeta[0]:
            failedFilters.append('AREAMAPQ')
        if 'AREASMS' in filters2Apply and not areaMeta[1]:
            failedFilters.append('AREASMS')

    ## 8. FILTER 8: Whether a metacluster has a SV_type assigned or not
    if 'IDENTITY' in filters2Apply: 
        if not identityFilter(metacluster, mode=mode):
            failedFilters.append('IDENTITY')
    
    ## 9. FILTER 9: Whether a metacluster is germline
    if 'GERMLINE' in filters2Apply:
        if not filter_germline_metaclusters(metacluster, confDict['minNormalReads']):
            failedFilters.append('GERMLINE')
    
    ## 10. FILTER 10: Reciprocal metaclusters range 
    if 'META-RANGE' in filters2Apply: 
        if not filter_clusterRange_reciprocalMeta(metacluster, 40):
            failedFilters.append('META-RANGE')
    
    ## 11. FILTER 11: Whether a cluster is on a L1PA element and has no pA support
    if 'ANNOTATION' in filters2Apply:
        if not metacluster.pA:
            if not filter_metacluster_subfamilyAnnot(metacluster, 'L1', 'L1PA', 0):
                failedFilters.append('ANNOTATION')
    
    ## 12. FILTER 12: Whether a cluster is on noisy region (calculated using the normal)
    if 'SVs-NORMAL' in filters2Apply:
        maxPercSVs = 20
        if not filter_percSVs_inNormal(metacluster, maxPercSVs, normalBam, confDict):
            failedFilters.append('SVs-NORMAL')
    
    ## 13. FILTER 13: Whether a transduction and the source element are in chr6 with no pA support
    if 'srcREF-TDs-ref6' in filters2Apply:
        if not filter_TDs_chr6_srcElement(metacluster):
            failedFilters.append('srcREF-TDs-ref6')
                
    return failedFilters

# def filter_metacluster_sonia(metacluster, filters2Apply, bam, normalBam, confDict):
#     '''
#     Apply selected filters to one metacluster.

#     Input:
#         1. metacluster: metacluster object
#         2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
#         3. bam  
#         4. normalBam
#         5. confDict

#     Output:
#         1. failedFilters -> list containing those filters that the metacluster doesn't pass.
#     '''        
#     failedFilters = []

#     ## 1. FILTER 1: Minimum number of reads per cluster
#     if 'MIN-NBREADS' in filters2Apply and not filter_min_nb_reads(metacluster, confDict['minReads'], confDict['minNormalReads']):
#         failedFilters.append('MIN-NBREADS')

#     ## 2. FILTER 2: Maximum number of reads per cluster
#     elif 'MAX-NBREADS' in filters2Apply and not filter_max_nb_reads(metacluster, confDict['maxClusterSize']):
#         failedFilters.append('MAX-NBREADS')
            
#     ## 3. FILTER 3: Area mapping quality
#     elif 'AREAMAPQ' in filters2Apply and not area(metacluster, confDict, bam)[0]:
#         failedFilters.append('AREAMAPQ')
    
#     ## 4. FILTER 4: Area clipping SMS  
#     elif 'AREASMS' in filters2Apply and not area(metacluster, confDict, bam)[1]:
#         failedFilters.append('AREASMS')

#     ## 5. FILTER 5: Whether a metacluster has a SV_type assigned or not
#     elif 'IDENTITY' in filters2Apply and not identityFilter(metacluster):
#         failedFilters.append('IDENTITY')
    
#     ## 6. FILTER 6: Whether a metacluster is germline
#     elif 'GERMLINE' in filters2Apply and not filter_germline_metaclusters(metacluster, confDict['minNormalReads']):
#         failedFilters.append('GERMLINE')
    
#     ## 7. FILTER 7: Reciprocal metaclusters range 
#     elif 'META-RANGE' in filters2Apply and not filter_clusterRange_reciprocalMeta(metacluster, 50):
#         failedFilters.append('META-RANGE')
    
#     ## 8. FILTER 8: Whether a cluster is on a L1PA element and has no pA support
#     elif 'ANNOTATION' in filters2Apply and not metacluster.pA and not filter_metacluster_subfamilyAnnot(metacluster, 'L1', 'L1PA', 100):
#         failedFilters.append('ANNOTATION')
    
#     ## 9. FILTER 9: Whether a cluster is on a L1PA element and has no pA support
#     elif 'SVs-NORMAL' in filters2Apply and not filter_percSVs_inNormal(metacluster, 30, normalBam, confDict):
#         failedFilters.append('SVs-NORMAL')

#     return failedFilters


def filter_min_nb_reads(cluster, minReads, minNormalReads):
    '''
    Filter cluster by comparing the number of supporting events with a minimum treshold

    Input:
        1. cluster: cluster object
        2. minReads: min number of events threshold
        3. minReads: min number of events threshold for normal sample

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    if cluster.events:
        
        ## 1. Compute number of events supporting the cluster 
        nbTotal, nbTumour, nbNormal = cluster.supportingReads()[0:3]

        ## 2. Compare the number of events supporting the cluster against the minimum required
        # 2.1 Paired mode:
        if nbTumour != None:

            if nbTumour >= minReads and nbNormal >= minNormalReads:
                cluster.mutOrigin = 'germline'
                PASS = True

            elif nbTumour >= minReads and not nbNormal >= minNormalReads:
                cluster.mutOrigin = 'somatic-tumour'
                PASS = True

            elif nbNormal >= minReads:
                cluster.mutOrigin = 'somatic-normal'
                PASS = True

            else:
                PASS = False

        # 2.1 Single mode:
        else:
            # If running in single mode (do no set mutation origin because it doesn't make sense)
            if nbTotal >= minReads:
                PASS = True
            else:
                PASS = False
                
    else:
        PASS = False    
    
    return PASS
    
    

def filter_max_nb_reads(metacluster, maxNbEvents):
    '''
    Filter metacluster by comparing the number of cluster supporting events with a maximum treshold

    Input:
        1. metacluster: metacluster object
        2. maxNbEvents: maximum number of events threshold
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 1. Compute number of events supporting the cluster 
    nbTotal = metacluster.nbEvents()[0]

    ## 2. Compare the number of events supporting the cluster against the maximum required
    if nbTotal <= maxNbEvents:
        PASS = True
    else:
        PASS = False
    
    return PASS

def filter_max_cv(metacluster, maxClusterCV):
    '''
    Filter metacluster by comparing its Coefficient of Variation with a maximum threshold.

    Input:
        1. metacluster: metacluster object
        2. maxClusterCV: maximum Coefficient of Variation threshold
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 1. Compute CV of the cluster 
    cv = metacluster.subclusters['INS'].cv_len()[1]

    ## 2. Compare the cluster CV against the maximum required
    if cv <= maxClusterCV:
        PASS = True
    else:
        PASS = False

    return PASS

def filter_perc_resolved(metacluster, minPercResolved):
    '''
    Filter metacluster by comparing the % of inserted sequence resolved with a minimum threshold

    Input:
        1. metacluster: metacluster object
        2. minPercResolved: minimum % resolved

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    if (metacluster.SV_features['PERC_RESOLVED'] >= minPercResolved):
        PASS = True

    else:
        PASS = False

    return PASS

def filter_SV_type(metacluster, targetSV):
    '''
    Filter metacluster by checking its SV type

    Input:
        1. metacluster: metacluster object
        2. targetSV: list containing list of target SV types
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    if metacluster.SV_type in targetSV:
        PASS = True 
        
    else:
        PASS = False

    return PASS

def identityFilter(cluster, mode ='SR'):
    '''
    Filter metacluster by checking its SV type

    Input:
        1. metacluster: metacluster object
        2. mode: SR or LR (Short Reads or Long Reads). Default='SR'
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    if cluster.identity != None:
        PASS = True 
        
    else:
        PASS = False

    if mode == 'LR':

        if 'IDENTITY' in cluster.SV_features.keys():
            if cluster.SV_features['IDENTITY']:
                PASS = True 
            else:
                PASS = False
        else:
            PASS = False  
    
    return PASS

def filter_suppl_MAPQ(clipping, minMAPQ):
    '''
    Filter out clipping cluster based on the average mapping quality of its supplementary alignments

    Input:
        1. clipping: clipping cluster object
        2. minMAPQ: minimum mapping quality

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't    
    '''

    ## Compute average mapping quality for supplementary alignments
    MAPQs = [event.mapQ for event in clipping.supplCluster.events]
    meanMAPQ = statistics.mean(MAPQs)

    ## Apply filter
    if meanMAPQ >= minMAPQ:
        PASS = True

    else:
        PASS = False
    
    return PASS

def filter_clipping_suppl_ref(clipping, targetRefs):
    '''
    Filter out clipping cluster if its supplementary cluster not located over a target reference

    Input:
        1. clipping: clipping cluster
        2. targetRefs: List of target references

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    ## Retrieve suppl. alignment reference 
    if clipping.supplCluster.ref in targetRefs:
        PASS = True 

    else:
        PASS = False

    return PASS

def filter_discordant_mate_ref(discordant, targetRefs):
    '''
    Filter out discordant cluster if its mate not located over a target reference

    Input:
        1. discordant: discordant cluster object
        2. targetRefs: List of target references

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    ## Retrieve mates referene 
    matesRef = discordant.events[0].mateRef

    if matesRef in targetRefs:
        PASS = True 
        
    else:
        PASS = False

    return PASS


def filter_clipping_suppl_position(clipping, ranges, buffer):
    '''
    Filter out clipping cluster if suppl. alignment within one of the provided regions
    
    Input:
        1. clipping: clipping cluster object
        2. ranges: Dictionary with reference ids as keys and the list of ranges on each reference as values
        3. buffer: Extend each range at their begin and end coordinate by a number of nucleotides == buffer length

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    ## Do not filter out clipping cluster if no input range on that particular reference 
    if clipping.supplCluster.ref not in ranges:
        PASS = True
        return PASS
        
    ## Assess overlap between mates interval and provided regions. 
    PASS = True

    for interval in ranges[clipping.supplCluster.ref]:
        rangeBeg, rangeEnd = interval

        # Add buffer
        rangeBeg = rangeBeg - buffer
        rangeEnd = rangeEnd + buffer

        # Assess overlap
        overlap, overlapLen = gRanges.overlap(clipping.supplCluster.beg, clipping.supplCluster.end, rangeBeg, rangeEnd)

        if overlap:
            PASS = False

    return PASS

def filter_discordant_mate_position(discordant, ranges, buffer):
    '''
    Filter out discordant cluster if mates align within one of the provided regions

    Input:
        1. discordant: discordant cluster object
        2. ranges: Dictionary with reference ids as keys and the list of ranges on each reference as values
        3. buffer: Extend each range at their begin and end coordinate by a number of nucleotides == buffer length

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## Retrieve mates position and interval
    matesRef = discordant.events[0].mateRef
    matesBeg, matesEnd = discordant.mates_start_interval()

    ## Do not filter out discordant cluster if no input range on that particular reference 
    if matesRef not in ranges:
        PASS = True
        return PASS
        
    ## Assess overlap between mates interval and provided regions. 
    PASS = True

    for interval in ranges[matesRef]:
        rangeBeg, rangeEnd = interval

        # Add buffer
        rangeBeg = rangeBeg - buffer
        rangeEnd = rangeEnd + buffer

        # Assess overlap
        overlap, overlapLen = gRanges.overlap(matesBeg, matesEnd, rangeBeg, rangeEnd)

        if overlap:
            PASS = False

    return PASS

def filter_discordant_mate_MAPQ(discordant, minMAPQ, bam, normalBam):
    '''
    Filter out discordant clusters based on average MAPQ for mate alignments

    Input:
        1. discordant: discordant cluster object
        2. minMAPQ: minimum average of mapping quality for mate alignments
        3. bam: path to bam file containing alignments for discordant cluster supporting reads and their mate
        4. normalBam: path to the matched normal bam file. If running in single mode, set to 'None' 

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    ## Open BAM files for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    if normalBam != None:
        bamFile_normal = pysam.AlignmentFile(normalBam, "rb")
        
    ## Define interval to search for mate alignment objects
    matesBeg, matesEnd = discordant.mates_start_interval()

    intervalRef = discordant.events[0].mateRef 
    intervalBeg = matesBeg - 500
    intervalEnd = matesEnd + 500

    ## Collect cluster supporting reads
    nbTotal, nbTumour, nbNormal, readIds, readIdsTumour, readIdsNormal = discordant.supportingReads()

    ## Compute average mapping quality for mates of cluster supporting reads
    # if running in single mode or running in paired but no reads from the normal are found in this interval
    if (normalBam == None or (normalBam != None and readIdsNormal == [])):
            
        avMAPQ = bamtools.average_MAPQ_reads_interval(intervalRef, intervalBeg, intervalEnd, readIds, bamFile)
                    
        if avMAPQ >= minMAPQ:
            PASS = True

        else:
            PASS = False
                
    # if running in paired mode and there is reads in the interval belonging to the normal bam
    else:
            
        avMAPQ = bamtools.average_MAPQ_reads_interval(intervalRef, intervalBeg, intervalEnd, readIdsTumour, bamFile)
        avMAPQ_normal = bamtools.average_MAPQ_reads_interval(intervalRef, intervalBeg, intervalEnd, readIdsNormal, bamFile_normal)
            
        # tumor and normal MAPQ average
        avMAPQ_pair = (avMAPQ * len(readIdsTumour) + avMAPQ_normal * len(readIdsNormal))/len(readIdsTumour + readIdsNormal)          
             
        if avMAPQ_pair >= minMAPQ:
            PASS = True

        else:
            PASS = False

    ## Close 
    bamFile.close()

    if normalBam != None:
        bamFile_normal.close()
        
    return PASS


def filter_germline(cluster, minNormal):
    '''
    Filter out those clusters formed by tumour and normal reads
    
    Input:
        1. cluster: cluster object
        2. minNormal: minimum number of reads supporting a SV in normal sample
    
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    count = 0
        
    for event in cluster.events:
        
        if event.sample == 'NORMAL':
            count += 1
                
    if count < minNormal:
        PASS = True
        
    else:
        PASS = False
        
    return PASS


def filter_germline_metaclusters(metacluster, minNormal):
    '''
    Filter out those clusters formed by tumour and normal reads
    
    Input:
        1. cluster: cluster object
        2. minNormal: minimum number of reads supporting a SV in normal sample
    
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    
    PASS = True
    
    nbNormal = metacluster.supportingReads()[2]
    
    if nbNormal >= minNormal:
        PASS = False
        
    return PASS



def filter_discordant_mate_unspecific(discordant, threshold, bam):
    '''
    Filter out discordant whose mates are located in regions captured unspecifically
    Insertion points where there is more than discordant reads are likely to be false positives
    Example of recurrent false positive filtered: chr6:29765954 (hg19)
    
    Input:
        1. discordant: discordant cluster instance
        2. threshold: ratio of nbDiscordant reads between nbTotal reads
        3. bam: path to bam file

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    
    # Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")
    
    nbProperPair = 0
    nbDiscordants, discordantList = discordant.nbReads()
        
    buffer = 200
    ref = discordant.events[0].mateRef
    beg, end = discordant.mates_start_interval()
        
    # Extract alignments
    iterator = bamFile.fetch(ref, beg - buffer, end + buffer)
        
    # Count properly paired reads in region
    properPairs = []
    for alignmentObj in iterator:
            
        if alignmentObj.is_proper_pair and not alignmentObj.is_supplementary:
                
                properPairs.append(alignmentObj.query_name)
                
    nbProperPair = len(set(properPairs))
        
    # if there are properly paired reads around insertion
    if nbProperPair > 0:
               
        if nbDiscordants/nbProperPair > threshold:        
            PASS = True

        else:  
            PASS = False

    else:
         
        PASS = True
            
    return PASS


# --------------- SHORT READS -----------------------
# HACER OTRA PARECIDA A LA QUE ESTABA PARA SHORT READS

## [SR CHANGE]
def applyFilters(clusters):
    '''
    Remove those clusters that fail in one or more filters
    '''
    newClusterDict = {}
    clusterTypes = clusters.eventTypes

    # NOTE: It would be better remove it from the list instead of create a new list
    # NOTE: It would be better get clusterType in another way, without performing two loops
    for clusterType in clusterTypes:
        for cluster in clusters.collect([clusterType]):
            newClusterDict[clusterType]=[]
            if False not in [value for value in cluster.filters.values()]:
                newClusterDict[clusterType].append(cluster)

    return newClusterDict


# [SR CHANGE]
def area(cluster,confDict,bam):
    '''
    Apply filters to cluster (SR bam) based on the characteristics of its region.

    Input:
        1. cluster: cluster object
        2. confDict
        3. bam
    Output:
        1. percMAPQFilter -> boolean: True if the cluster pass the filter, False if it doesn't
        2. percSMSReadsFilter -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## Set region coordinates
    if cluster.beg == cluster.end:
        if cluster.beg > 100:
            binBeg = cluster.beg - 100
        else:
            binBeg = cluster.beg
        binEnd = cluster.beg + 100
    else:
        if cluster.beg > 50:
            binBeg = cluster.beg - 50
        else:
            binBeg = cluster.beg
        binEnd = cluster.end + 50

    ref = cluster.ref

    ## Extract filter parameters from config dict
    minReadsRegionMQ = confDict['minReadsRegionMQ']
    maxRegionlowMQ = confDict['maxRegionlowMQ']
    maxRegionSMS = confDict['maxRegionSMS']

    # Set counts to 0
    lowMAPQ = 0
    SMSReads = 0
    nbReads = 0

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, binBeg, binEnd)

    for alignmentObj in iterator:
        
        if alignmentObj.cigartuples != None:
        
            # Check if aligment pass minimum mapq for reads within the cluster region
            passMAPQ = areaMAPQ(alignmentObj, minReadsRegionMQ)

            # If it doesnt pass, add 1 to the counts of low mapping quality reads within the cluster region
            if passMAPQ == False:
                lowMAPQ += 1

            # Check if aligment is mapped this way: Soft Match Soft (SMS)
            SMSRead = areaSMS(alignmentObj)
            
            # If it is mapped SMS, add 1 to the counts of SMS reads within the cluster region
            if SMSRead == True:
                SMSReads += 1
            # Count total number of reads in the region
            nbReads += 1
    
    ## Calculate percentages
    percMAPQ = stats.fraction(lowMAPQ, nbReads)
    percSMSReads = stats.fraction(SMSReads, nbReads)

    ## If the percentage of low MQ reads is lower than the threshold pass the filter.
    if percMAPQ != None:
        if percMAPQ < float(maxRegionlowMQ):
            percMAPQFilter = True
        else:
            percMAPQFilter = False
    else:
        percMAPQFilter = False

    ## If the percentage of SMS reads is lower than the threshold pass the filter.
    if percSMSReads != None:
        if percSMSReads < float(maxRegionSMS):
            percSMSReadsFilter = True
        else:
            percSMSReadsFilter = False
    else:
        percSMSReadsFilter = False
    
    return percMAPQFilter, percSMSReadsFilter

# [SR CHANGE]
def areaMAPQ(alignmentObj, minReadsRegionMQ):
    '''
    Check if the MAPQ of a read pass the minReadsRegionMQ threshold

    Input:
        1. alignmentObj
        2. minReadsRegionMQ
    Output:
        1. percMAPQFilter -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    MAPQ = int(alignmentObj.mapping_quality)

    if MAPQ > int(minReadsRegionMQ):
        passMAPQ = True
    else:
        passMAPQ = False

    return passMAPQ

# [SR CHANGE]
def areaSMS(alignmentObj):
    '''
    Check if aligment is mapped this way: Soft Match Soft (SMS)

    Input:
        1. alignmentObj
        2. maxRegionSMS
    Output:
        1. percSMSReadsFilter -> boolean: True if the cluster pass the filter, False if it doesnt
    '''

    SMSRead = False

    # Select first and last operation from cigar to search for clipping
    firstOperation, firstOperationLen = alignmentObj.cigartuples[0]
    lastOperation, lastOperationLen = alignmentObj.cigartuples[-1]
    ## Clipping >= X bp at the left
    #  Note: soft (Operation=4) or hard clipped (Operation=5)     
    if ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation == 4) or (lastOperation == 5)):
        SMSRead = True
        
    #if ((lastOperation == 4) or (lastOperation == 5)) and ((firstOperation != 4) and (firstOperation != 5)):
        #SMSRead = True

    return SMSRead


def filter_highDup_clusters(cluster, maxDupPerc):
    '''
    Filter out those clusters formed by more than a percentage of duplicates
    
    Input:
        1. cluster: list of discordant clusters formed by DISCORDANT events
        2. maxDupPerc: Maximun % of duplicates allowed (not included)
    
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    
    dupPerc = cluster.dupPercentage()
        
    if dupPerc < maxDupPerc:
            
        PASS = True

    else:
        PASS = False
                
    return PASS



def filter_clusterRange_discordant(cluster):
    '''
    Filter out those discordant clusters in which all reads are piled up.
    This filter is only applied to clusters formed by more than a discordant alignment
    
    ----------****>                       ---------******>
       -------*******>                    ---------******>
     ---------*****>                      ---------******>
         -----*********>                  ---------******>
    |         |                          |         |   
    beg      end                         beg      end
    ----------       clusterRange         ----------
         -----      min(alignRanges)      ----------
       True             PASS                 False      
    
    Input:
        1. cluster: cluster formed by DISCORDANT events
    
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    
    PASS = True
        
    # if there is more than a discordant alignment
    if len(cluster.events) > 1:
               
        # define cluster range
        clusterRange = cluster.end - cluster.beg
        
        # define minimum alignment range of all reads supporting the cluster
        readRanges = []
        
        for event in cluster.events:
            beg, end = event.readCoordinates()
            readRange = abs(abs(end) - abs(beg))
            readRanges.append(readRange)
        
        alignRange = min(readRanges)
        
        # if the cluster range is smaller or equal to the minimum alignment range
        if (clusterRange <= alignRange):
            
            # discard the cluster
            PASS = False
    
    return PASS



def filter_clusterRange_clipping(cluster):
    '''
    Filter out those clipping clusters in which all clippings have the same 
    coordinates relative to the read. This filter is only applied to clusters 
    formed by more than a clipping alignment
    
    ----------****>                       ---------******>
             55   75                              40    75
       -------*******>                    ---------******>
              38    75                            40    75
     ---------*****>                      ---------******>
              61   75                             40    75
                          clusterCoord 
    [55, 75, 38, 75, 61, 75]        [40, 75, 40, 75, 40, 75]
          True              PASS                False      
    
    Input:
        1. cluster: cluster formed by CLIPPING events
    
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    
    PASS = True
       
    # if there is more than a clipping
    if len(cluster.events) > 1:
               
        # with clipping events beg and end coordinates are the same. However, readCoordinates()
        # returns the aligment coordinates relative to the read. If more than 2 coordinates, 
        # it is not a piled up cluster of clippings        
        clusterCoord = [event.readCoordinates() for event in cluster.events]
        clusterCoord_flatList = [item for sublist in clusterCoord for item in sublist]
        n_clusterCoord = len(set(clusterCoord_flatList))
              
        if(n_clusterCoord <= 2):
            PASS = False
        
    return PASS


def filter_clusterRange_reciprocalMeta(metacluster, maxPercOverlap):
    '''
    Filter metacluster reciprocal clusters if plus clusters and minus clusters overlap more than buffer
    Filter out:
         ------->
    <-------
     ------->
      <-------
             <-------
   end         beg
   
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    PASS = True
    
    # if reciprocal cluster
    if metacluster.orientation == "RECIPROCAL":
        
        # create subclusters
        subclusters = metacluster.create_subclusters()
        
        # if discordant clusters formed
        if 'MINUS-DISCORDANT' in subclusters.keys() and 'PLUS-DISCORDANT' in subclusters.keys():
            
            # determine subclusters coordinates
            minus_beg, minus_end = subclusters['MINUS-DISCORDANT'].beg, subclusters['MINUS-DISCORDANT'].end
            plus_beg, plus_end = subclusters['PLUS-DISCORDANT'].beg, subclusters['PLUS-DISCORDANT'].end
            
            # determine overlap between minus and plus clusters
            overlap, overlapLen = gRanges.rcplOverlap(minus_beg, minus_end, plus_beg, plus_end, 1)

            # if overlap
            if overlap:
                
                # calculate ratio of overlapLen / total metacluster amplitude 
                lenCluster = metacluster.end - metacluster.beg
                percOverlap = overlapLen/lenCluster * 100
                
                if percOverlap > maxPercOverlap:
                    PASS = False
        
    return PASS


def filter_metacluster_subfamilyAnnot(metacluster, family, subtype, minDist):
    '''
    Filter metacluster if located at less than minDist in specific ME subtypes
    
    Input:
    1. metacluster
    2. family: Annotation family. Ex: 'L1'
    3. subtype: Annotation subfamily. Ex: 'L1PA'
    4. minDist: minimum distance from metacluster to annotated repeat. Ex: 150
    
    Output: 
    1. PASS: Boolean indicating whether it should be filtered out or not
    '''
    PASS = True
    
    if metacluster.repeatAnnot:
        
        repeats = metacluster.repeatAnnot
        families = [repeat['family'] for repeat in repeats]
        
        if family in families:
            
            subfamilies = [repeat['subfamily'] for repeat in repeats] 
            distances = [repeat['distance'] for repeat in repeats]
            
            for idx, subfamily in enumerate(subfamilies):
                
                if subtype in subfamily:
                    
                    if distances[idx] <= minDist:
                        
                        PASS = False
                        break
        
    return(PASS)


def filter_percSVs_inNormal(metacluster, maxPercSVs, normalBam, confDict):
    '''
    Filter metacluster if region in normal has a percentage of discordants greater than maxPercDisc
    
    Input:
        1. metacluster
        2. maxPercSVs: maximum percentage of discordants in region 
        3. normalBam
        4. confDict
    
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    PASS = True
    
    # define beg and end coordinates
    beg = metacluster.refLeftBkp if metacluster.refLeftBkp is not None else metacluster.beg
    end = metacluster.refRightBkp if metacluster.refRightBkp is not None else metacluster.end
    
    # count nb of reads pointing to SVs
    confDict['minMAPQ'] = 1
    confDict['targetEvents'] = ['DISCORDANT', 'CLIPPING']
    eventsDict = bamtools.collectSV(metacluster.ref, beg, end, normalBam, confDict, None, supplementary = True)
    
    nbSVs = 0
    if 'DISCORDANT' in eventsDict.keys():
        nbSVs += len(eventsDict['DISCORDANT'])
    if 'LEFT-CLIPPING' in eventsDict.keys():
        nbSVs += len(eventsDict['LEFT-CLIPPING'])
    if 'RIGHT-CLIPPING' in eventsDict.keys():
        nbSVs += len(eventsDict['RIGHT-CLIPPING'])
    
    # count total number of reads in region 
    bamFile = pysam.AlignmentFile(normalBam, "rb")
    nbReads = bamFile.count(metacluster.ref, beg, end)
    
    if nbReads > 0:
        percSVs = nbSVs/nbReads * 100
        
        if percSVs > maxPercSVs:
            PASS = False
            
    return(PASS)
     

def filter_TDs_chr6_srcElement(metacluster):
    '''
    Filter metacluster if there is a transduction in chr 6 and the source element is located in chr 6 as well.
    It is a source of FPs. Metalusters will only be discarded if no polyA support is found.
    
    Input:
    1. metacluster
    
    Output: 
    1. PASS: Boolean indicating whether it should be filtered out or not
    '''
    PASS = True
    
    # If it's a transduction in chr6
    if metacluster.ref == '6' and metacluster.identity in ['TD1', 'TD2', 'partnered', 'orphan']:
        
        # Extract src_id ref
        src_id = metacluster.src_id.replace('p', 'q')
        src_ref = src_id.split('q')[0]
        
        print("filter_TDs_sameRef_srcElement")
        print("src_ref, metacluster.ref, metacluster.pA")
        print(str(src_ref))
        print(str(metacluster.ref))
        print(str(metacluster.pA))
        
        # If scr_ref is equal to transduction ref and there is no pA support
        if metacluster.ref == src_ref and not metacluster.pA:
            
            PASS = False
        
    return(PASS)