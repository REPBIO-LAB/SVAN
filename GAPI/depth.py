'''
Module 'callers' - Contains classes and functions for calling variants from next generation sequencing data
'''

## DEPENDENCIES ##
# External
import multiprocessing as mp
import os
import time
import pysam
import numpy as np

# Internal
from GAPI import bamtools 

## FUNCTIONS ##


## CLASSES ##
class segment():
    '''
    Genomic segment
    '''
    def __init__(self, ref, beg, end, bam, sampleId):
        '''
        Initialize object instance
        
        Input:
            1. ref: reference
            2. beg: begin position
            3. end: end position
            4. bam: path to bam file
            5. sampleId: sample identifier
        '''
        ## General:
        self.ref = ref
        self.beg = beg
        self.end = end
        self.bam = bam
        self.sampleId = sampleId
        
        ## Read depth:
        self.TDP = 0 # Total 
        self.WDP = 0 # Watson (+)
        self.CDP = 0 # Crick (-)

        ## Read depth ratios/fold changes:
        self.TDP_fc = None
        self.WDP_fc = None
        self.CDP_fc = None

    def coordId(self):
        '''
        Return segment coordinates identifier as ref:beg-end
        '''
        return self.ref + ':' + str(self.beg) + '-' + str(self.end)

    def read_depth(self, bam, minMAPQ, filterDup):
        '''
        Compute read depth on the segment. 

        Input:
            1. bam: path to bam file
            2. minMAPQ: minimum read mapping quality
            3. filterDup: filter read duplicates (True) or not (False)

        Output: 
            Update TDP, WDP, CDP counters
        '''
        ## Open BAM file
        bamFile = pysam.AlignmentFile(bam, "rb")
    
        ## Extract alignments
        iterator = bamFile.fetch(self.ref, self.beg, self.end)
    
        # For each read alignment
        for alignmentObj in iterator:
        
            ### 1. Alignment filtering:
            ## a) No query sequence available
            if alignmentObj.query_sequence == None:
                continue

            ## b) Unmapped read 
            if alignmentObj.is_unmapped:
                continue

            ## c) Minimum mapping quality
            MAPQ = int(alignmentObj.mapping_quality) 

            if (MAPQ < minMAPQ):
                continue                

            ## d) Duplicates filtering enabled and duplicate alignment
            if filterDup and alignmentObj.is_duplicate:
                continue
        
            ## e) Supplementary alignment 
            if alignmentObj.is_supplementary:
                continue

            ## f) Secondary alignment
            if alignmentObj.is_secondary:
                continue
        
            ### 2. Update counters
            self.TDP += 1

            # a) Watson (+)
            if not alignmentObj.is_reverse:
                self.WDP += 1
        
            # b) Crick (-)
            else:
                self.CDP += 1        
        
        ## Close bam
        bamFile.close()

class consensus_segment():
    '''
    Consensus genomic segment
    '''
    def __init__(self, ref, beg, end, segments):
        '''
        Initialize object instance
        
        Input:
            1. ref: reference
            2. beg: begin position
            3. end: end position
            4. segments: list of segment instances spanning the same genomic position
        '''
        self.ref = ref
        self.beg = beg
        self.end = end
        self.segments = segments

        ## Read depth attributes
        self.TDP = None # Total 
        self.WDP = None # Watson (+)
        self.CDP = None # Crick (-)
        self.TDP_sd = None # Total 
        self.WDP_sd = None # Watson (+)
        self.CDP_sd = None # Crick (-)

    def coordId(self):
        '''
        Return segment coordinates identifier as ref:beg-end
        '''
        return self.ref + ':' + str(self.beg) + '-' + str(self.end)

    def consensus_read_depth(self):
        '''
        Compute compute consensus read depth by taking into account the depth across all segments

        Output: 
            Define TDP, WDP, CDP and their standard deviations (TDP_sd, WDP_sd, CDP_sd)
        '''
        ## 1. Total read depth
        TDP_list = [segment.TDP for segment in self.segments]
        self.TDP = np.mean(TDP_list)
        self.TDP_sd = np.std(TDP_list)

        ## 2. Watson read depth
        WDP_list = [segment.WDP for segment in self.segments]
        self.WDP = np.mean(WDP_list)
        self.WDP_sd = np.std(WDP_list)

        ## 3. Crick read depth
        CDP_list = [segment.CDP for segment in self.segments]
        self.CDP = np.mean(CDP_list)
        self.CDP_sd = np.std(CDP_list)

class read_depth_caller():
    '''
    Read depth caller 
    '''
    def __init__(self, bam, sampleId, confDict):
        '''
        Initialize object instance

        Input:
            1. bam: input bam file
            2. sampleId: sample identifier
            3. confDict: configuration dictionary with the following key value pairs:
                - binSize: segments size 
                - targetRefs: list with target references
                - processes: number of processes for parallelization
                - minMAPQ: minimum read mapping quality
                - filterDup: filter read duplicates (True) or not (False)
        '''
        self.bam = bam
        self.sampleId = sampleId
        self.confDict = confDict

        ## Compute reference lengths
        self.refLengths = bamtools.get_ref_lengths(self.bam)

    def read_depth_wg(self):
        '''
        Compute read depth genome wide across non-overlapping bins or for a predefined set of intervals
        '''
        ### 1. Define genomic bins for read depth computation ##
        bins = bamtools.binning(None, self.bam, self.confDict['binSize'], self.confDict['targetRefs'])

        ## Select the first and last 10 bins for testing:
        #bins = bins[:11] + bins[-9:]
        
        ### 2. Calculate read depth per genomic bin ##
        segments = [self.read_depth_segment(ref, beg, end) for ref, beg, end in bins]

        return segments

    def read_depth_wg_mp(self):
        '''
        Compute read depth genome wide across non-overlapping bins or for a predefined set of intervals. 
        Computation is speed up via multi-processing
        '''
        ### 1. Define genomic bins for read depth computation ##
        bins = bamtools.binning(None, self.bam, self.confDict['binSize'], self.confDict['targetRefs'])

        ## Select the first and last 10 bins for testing:
        bins = bins[:11] + bins[-9:]

        ### 2. Calculate read depth per genomic bin ##
        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        segments = pool.starmap(self.read_depth_segment, bins)
        pool.close()
        pool.join()

        return segments

    def read_depth_segment(self, ref, beg, end):
        '''
        Compute read depth for a genomic segment/bin
        '''
        ## Create segment object
        segmentObj = segment(ref, beg, end, self.bam, self.sampleId)

        ## Compute read depth
        segmentObj.read_depth(self.bam, self.confDict['minMAPQ'], self.confDict['filterDup'])

        ## Return segment
        return segmentObj