'''
Module 'singleCell' - Contains classes for with single cell data (SVs, for now...)
'''

# External
import random    
import copy
import itertools
from scipy import stats
import statistics
import math 
from statsmodels.stats.multitest import multipletests
import pandas as pd
import operator
import collections
import numpy as np

# Internal 
import gRanges
import structures
import clustering
import variants
import output

class population():
    '''
    Class representing a population of single cells
    '''
    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.cells = {}

    def create_cells(self, cellIds, refLen):
        '''
        Create set of cells 

        Input:
            1. cellIds: list of cell identifiers
            2. refLen: dictionary with a key per reference/chromosome and its length as value

        Output:
            Update 'cells' argument with a dictionary containing cellIds as keys and cell instances as values
        '''
        ## For each cell
        for cellId in cellIds:
            
            ## Initialize cell
            cellInstance = cell(cellId, refLen)

            ## Add to the dict
            self.cells[cellId] = cellInstance

    def nbCells(self):
        '''
        Return the number of cells in the population
        '''
        return len(self.cells)

    def SVs2cells(self, SVs, refLen):
        '''
        Organize SV calls by cell

        Input:
            1. SVs: nested dictionary containing SVs identified with MosaiCatcher. One key per cellId and reference with the corresponding list of 'SV_call' objects as values
            2. refLen: dictionary with a key per reference/chromosome and its length as value

        Output:
            Update 'cells' argument with a dictionary containing cellIds as keys and cell instances as values
        '''

        ## 1. Initialize cells
        cellIds = list(SVs.keys())
        self.create_cells(cellIds, refLen)
            
        ## 2. Add SVs to the corresponding cell and chromosome
        # For each cell 
        for cellId, cell in self.cells.items():

            ## Add SVs
            cell.add_SV(SVs[cellId])

    def SV2segments(self):
        '''
        Aggregate all SVs identified across multiple cells into segments. Each segment will contain all identified
        SVs with begin and end (== segment). Each SV within a segment will derive from a single cell. 

        Output:
            1. segmentsOut: dictionary containing chromosome id as key and list of segment instances as values
        '''        
        ## 1. Agregate SVs with same begin and end coordinates into segments
        segments = {}

        ## For each cell
        for cellId, cell in self.cells.items():

            ## For each chromosome
            for chromId, chrom in cell.chromosomes.items():

                if chromId not in segments:
                    segments[chromId] = {}

                ## For each SV
                for SV in chrom.SVs:

                    segmentId = '_'.join([chromId, str(SV.beg), str(SV.end)])

                    ## Initialize segment
                    if segmentId not in segments[chromId]:
                        segments[chromId][segmentId] = variants.SV_segment(chromId, SV.beg, SV.end)
                        
                    ## Add SV to segment
                    segments[chromId][segmentId].add([SV])

        ## 2. Rearrange segments into a simple dictionary with chrom ids as keys
        segmentsOut = {}

        for chromId in segments:
            segmentsOut[chromId] = list(segments[chromId].values())

        return segmentsOut

    def correct_oversegmentation(self, refLen):
        '''
        Correct oversegmentation of SV calls by merging adjacent SVs with enough evidences supporting
        that they actually correspond to the same SV. Merging done at single-cell level
        
        Input:
            1. refLen: dictionary with a key per reference/chromosome and its length as value

        Output: Updated list of SVs per cell and chromosome after doing SV merging
        '''
        ## Merge adjacent SVs
        # For each cell
        for cellId, cell in self.cells.items():

            ## For each chromosome
            for chromId, chrom in cell.chromosomes.items():

                ## Do merging
                chrom.merge_SV()

    def consensus_SVs(self, refLen):
        '''
        For each SV compute consensus SV features taking into all the SVs identified across multiple cells 
        spanning the same segment

        Input:
            1. refLen: dictionary with a key per reference/chromosome and its length as value

        Output:
            For each SVs within a cell and chromosome add the following consensus attributes:
                - cAF: consensus allele frequency
                - cHaplo: list of consensus haplotypes 
                - cSV_type: list of consensus sV types
            segments
        '''
        ## 1. Aggregate SVs into segments
        segments = self.SV2segments()

        ## 2. Compute consensus segment features
        for chrom in segments:

            for segment in segments[chrom]:

                ## Consensus allele frequency (AF)
                segment.consensus_AF()

                ## Consensus SV type
                segment.consensus_SV_type()

                ## Consensus haplotype
                segment.consensus_haplotype() 
        
        return segments

     
    def search4clusteredSVs(self):
        '''
        For each cell in the population search for clusters of SVs

        Ouput:
        '''
        ## 1. Retrieve genomic windows per cell
        windows = {}

        for cellId, cell in self.cells.items():
            
            windows[cellId] = cell.chrom2windows_SVs(10000000, 3)

        ## 2. Organize all windows into a single dict (windowId as key)
        windows = structures.nestedDict2list(windows)
        windows = {window.id: window for window in windows}

        ## 3. Search for windows with SV clustering signatures
        ## Collect window features
        data = [[window.id, window.ref, window.nbSV(), window.median_SV_distance()] for window in windows.values()]

        ## Add header
        data.insert(0, ['name', 'chrom', 'nbSV', 'distance'])

        ## Create dataframe
        df = pd.DataFrame(data[1:], columns=data[0]).set_index('name')
        
        ## Cluster windows
        nbClusters, data = clustering.DBSCAN_clustering(df, 'nbSV', 'distance', 1, 4)
        
        ## Compute statistics per cluster
        header = ['clusterId', 'size', 'nbSV', 'dist']
        clusters = [header]
        
        for clusterId in range(0, nbClusters):

            ## Retrieve windows
            cluster = data[data['cluster'] == clusterId]

            ## Compute cluster stats
            size = len(cluster.index)  # == nb. windows
            medianNbSV = statistics.median(cluster['nbSV'])
            medianDist = statistics.median(cluster['distance'])

            ## Add stats
            clusters.append([clusterId, size, medianNbSV, medianDist])
            
        ## Convert to dataframe
        clusters = pd.DataFrame(clusters[1:], columns=clusters[0]).set_index('clusterId')

        ## Filter clusters based on two criteria:
        # - Medium distance between SVs across windows (threshold)
        # - Medium number of SVs across windows (threshold)
        filteredClusters = clusters[(clusters['dist'] <= 0) & (clusters['nbSV'] >= 4)]

        ## Pick windows composing each cluster
        for clusterId in filteredClusters.index.values:

            cluster = data[data['cluster'] == clusterId]        
            windowIds = cluster.index.values
            windowsInCluster = [windows[windowId] for windowId in windowIds]

    def search4complex(self, outDir):
        '''
        For each cell in the population search for chromosomes with complex rearrangements
        '''        
        ## 1. Assess each cromosome within every cell 
        BFBs = {}
        CTSs = {}

        # For each cell 
        for cellId, cell in self.cells.items():
            BFBs[cellId], CTSs[cellId] = cell.search4complex()

        ## 2. Cluster complex SVs identified across cells
        ## BFB
        BFB_clusters = variants.cluster_BFB(BFBs)

        ## Chromothripsis (CTS)
        CTS_clusters = variants.cluster_CTS(CTSs)

        ## 3. Filter out BFB based on haplotype consistency
        BFB_clusters = variants.haploFilter_BFB(BFB_clusters)

        ## 4. Report complex SVs
        ## Report BFBs single cell level calls
        output.report_BFB_per_cell(BFB_clusters, outDir)

        ## Report BFBs population level calls
        output.report_BFB(BFB_clusters, self.nbCells(), outDir)

        ## Report chromothripsis single cell level calls
        #output.report_CTS_per_cell(CTS_clusters, self.nbCells(), outDir)

        ## Report chromothripsis population level calls
        output.report_CTS(CTS_clusters, self.nbCells(), outDir)
        
class cell():
    '''
    Class representing a single cells
    '''
    def __init__(self, cellId, chrLen):
        '''
        Initialize class instance
        '''
        self.id = cellId

        ## Create chromosomes
        self.chromosomes = self.create_chromosomes(chrLen)

    def chrom_lengths(self):
        '''
        Return dictionary with chromosome ids as keys and lengths as values
        '''        
        chromLen = dict([(chrom.id, chrom.len) for chrom in self.chromosomes.values()])
        return chromLen
            
    def create_chromosomes(self, chrLen):
        '''
        Create chromosome instances for the cell

        Input:
            1. chrLen: dictionary with a key per chromosome and its length as value
        
        Output:
            1. chromosomes: Dictionary containing chromosome ids as keys and the corresponding chromosome instance as valuew
        '''
        chromosomes = {}

        ## For each chromosome
        for chrId, length in chrLen.items():

            # Create chromosome
            chrom = chromosome(chrId, length)
            chromosomes[chrom.id] = chrom  

            # Add cell id to the chromosome
            chromosomes[chrom.id].cellId = self.id
        
        return chromosomes

    def add_SV(self, SVs):
        '''
        Add SVs to cell's chromosomes

        Input:
            1. SVs: Dictionary with chromosome ids as keys and the list of SVs in each chromosome as values
        
        Output:
            Updated 'SVs' attribute for each chromosome in the cell
        '''    
        ## For each chromosome
        for chrId in SVs:

            ## Add SVs
            self.chromosomes[chrId].add_SV(SVs[chrId])      

    def search4complex(self):
        '''
        For each chromosome in the cell assess if it has complex SV signatures. 

        Several signatures are assessed:
            1. SV clustering. SVs not randomly distributed along the chromosome. Lower SV distance than expected by chance
            2. SV enrichment. Chromosome enriched in SVs, higher number of SVs than expected by chance
            3. Haplotype bias. SVs are biased towards a certain haplotype/s

        Output: Per chromosome update the following attributes:
            1. SV clustering signature. 'clustering' and 'c_pvalue'
            2. SV enrichment. 'enrichment' and 'e_pvalue'
            3. Haplotype bias. 'haploBias' and 'h_pvalue' 
        '''
        ## 1. Search for BFB
        BFBs = self.search4BFB()

        ## 2. Search for chromothripsis
        CTSs = self.search4CTS()

        return BFBs, CTSs
    

    def SV_clustering(self, minNbSV, nbRandom):
        '''
        Assess SV clustering per chromosome in the cell
        
        Input:
            1. minNbSV: minimum number of SVs in the chrom to assess clustering
            2. nbRandom: number of randomizations for computing random distribution

        Output: 
            Update three attributes in each chromosome
            1. clustering: boolean. Cluster with SV clustering (True) or not (False) 
            2. c_score: SV clustering score
            3. c_pvalue: p-value for clustering
        '''
        ## For each chromosome
        for chrom in self.chromosomes.values():

            chrom.SV_clustering(5, 10)   

    def SV_enrichment(self):
        '''
        Assess SV enrichment per chromosome in the cell

        Output: 
            Update three attributes in each chromosome
            1. enrichment: boolean. Cluster with SV clustering (True) or not (False) 
            2. e_score: SV enrichment score (+: enrichment; -: depletion)
            3. e_pvalue: p-value for enrichment  
        '''
        ## 1. Compute observed distribution (== number of SVs per chromosome)
        observed = dict([(chrom.id, chrom.nbSV()) for chrom in self.chromosomes.values()])

        ## 2. Compute random distribution 
        ## 2.1 Map chromosomes into 0-1 space
        chromLen = self.chrom_lengths()
        mappedRanges, ranges2chrom = gRanges.mapChrom(chromLen)

        ## 2.2 Randomly assign SVs to chromosomes
        ## Initialize counts
        expected = dict([(chrom.id, 0) for chrom in self.chromosomes.values()])

        ## Total number of SVs
        nbSV = sum([chrom.nbSV() for chrom in self.chromosomes.values()])

        ## Perform randomizations
        nbRandom = 100
        expected = []

        # Perform randomizations
        for i in range(nbRandom):
            randomization = self.assignSV2chrom(nbSV, mappedRanges, ranges2chrom)
            expected.append(randomization)

        # Group randomization counts by chromosome
        expected = structures.merge_dictionaries(expected)
        
        ## 3. Compare observed versus expected distribution per chromosome
        # For each chromosome
        for chrId, chrom in self.chromosomes.items():

            ## Fit expected data to normal distribution
            mean, std = stats.norm.fit(expected[chrId])

            ## Obtain enrichment/depletion significance
            chrom.e_pvalue = stats.norm.pdf(observed[chrId], loc=mean, scale=std)

            ## Compute enrichment score 
            pseudocount = 0.0000001
            enrichment_score = float(observed[chrId] + pseudocount) / mean + pseudocount 
            chrom.e_score = math.log(enrichment_score, 2)

            ## Assess for significance
            chrom.enrichment  = True if chrom.e_pvalue < 0.05 else False
            
    def haplo_bias(self):
        '''
        Assess haplotype bias per chromosome in the cell

        Output: 
            Update three attributes in each chromosome:
            1. haploBias: boolean. Cluster with haplotype bias (True) or not (False) 
            2. h_score: haplotype bias score (fisher exact test odd ratio)
            3. h_pvalue: p-value for haplotype bias       
        '''
        ## For each chromosome
        for chrom in self.chromosomes.values():
            chrom.haplo_bias()  

    def search4telDel(self):
        '''
        Search for telomeric deletions per chromosome in the cell
        '''
        ## For each chromosome
        for chrom in self.chromosomes.values():
            
            chrom.telomeric_del()  

    def cluster_SV(self):
        '''
        Cluster SVs per chromosome in the cell
        '''
        ## For each chromosome
        for chrom in self.chromosomes.values():
            chrom.cluster_SV()    

    def search4BFB(self):
        '''
        Search for BFB in each chromosome
        '''
        BFBs = {}

        ## For each chromosome
        for chrom in self.chromosomes.values():
            BFBs[chrom.id] = [chrom.search4BFB('p'), chrom.search4BFB('q')]
        
        return BFBs

    def search4CTS(self):
        '''
        Search for chromothripsis in each chromosome

        Output:
            1. chromothripsis: dictionary containing chromothriptic chromosomes. 
                               Keys are chromosome identifiers while values are
                               chromosome instances
        '''

        ### 1. Assess three metrics per chromosome
        ## 1.1 SV clustering
        self.SV_clustering(5, 10)

        ## 1.2 SV enrichment
        self.SV_enrichment()

        ## 1.3 Haplotype bias
        self.haplo_bias()

        ## 2. Assess for chromothripsis per chromosome 
        CTSs = {}

        for chrom in self.chromosomes.values():
            if chrom.is_chromothriptic():
                CTSs[chrom.id] = chrom
        
        return CTSs

    def assignSV2chrom(self, nbSV, ranges, ranges2chrom):
        '''
        Randomly asign a set of SVs to chromosomes. 

        Input:
            1. nbSV: number of SVs to be randomly assigned
            2. ranges: chromosomal ranges in 0-1 space (needed for making the probability dependent of the
                       length of the chromosome)
            3. ranges2chrom: mapping between chromosomal ranges and chromosome ids

        Output:
            1. expected: dictionary containing the number of randomly asigned SVs (value) to each chromosome (key)
        '''
        # Initialize counts
        expected = dict([(chrom.id, 0) for chrom in self.chromosomes.values()])

        # Randomly asign a SV to one chromosome
        for i in range(nbSV):

            ## Randomly pick one range
            randomNb = random.uniform(0, 1)
            index = gRanges.overlap_multiranges(randomNb, ranges)[1][0]
            targetRange = ranges[index]

            ## Range to chromosome mapping
            rangeId = '-'.join(map(str, targetRange))
            targetChr = ranges2chrom[rangeId]

            ## Add SV to cromosome by incrementing counter
            expected[targetChr] += 1

        return expected

    def chrom2windows_SVs(self, windowSize, minNbSV):
        '''
        For each chromosome retrieve genomic windows of length == windowSize. Windows will contain 
        the list of SVs falling within.  

        Input:
            1. windowSize: window size
            2. minNbSV: minimum number of SVs in a window to be reported

        Output: 
            1. windows: dictionary with chromosomes as keys and the list of windows per chromosome as value
        '''
        windows = {}

        ## Create windows per chromosome        
        for chrom in self.chromosomes.values():
            windows[chrom.id] = chrom.create_windows_SVs(windowSize, minNbSV)   

        return windows

    def merge_SV(self):
        '''
        Correct oversegmentation per chromosome by doing SV merging
        '''
        ## Merge SVs for each chromosome        
        for chrom in self.chromosomes.values():
            chrom.merge_SV()           

class chromosome():
    '''
    Class representing a chromosome
    '''
    def __init__(self, chrId, length):
        '''
        Initialize class instance
        '''
        self.id = chrId
        self.beg = 0
        self.end = length
        self.len = length
        self.SVs = []
        self.SV_clusters = []
        self.cellId = None

        ## Clustering parameters
        self.clustering = None 
        self.c_pvalue = None

    def add_SV(self, SVs):
        '''
        Add a set of SV_call events to the chromosome

        Input:
            1. SVs: list of SV_call objects

        Output: 
            Update SVs attribute with new SV_call events
        '''
        self.SVs = self.SVs + SVs

    def sort_SV(self):
        '''
        Sort SVs along the chromosome by their starting position (increasing)
        '''            
        self.SVs.sort(key=lambda SV: SV.beg, reverse=False)

    def nbSV(self):
        '''
        Compute the number of SVs in the chromosome

        Output:
            1. nbSVs: integer with the number of SVs
        '''
        return len(self.SVs)

    def SV_distance(self):
        '''
        Compute the distances between adjacent SVs along the chromosome

        Output:
            1. distances: list containing the distance between the begin and end position of adjacent SVs
        '''
        ## 1. Make sure that SVs are sorted by their start position
        self.sort_SV()

        ## 2. Compute distances between adjacent SVs
        distances = []

        ## For each SV and adjacent SV pair
        for SV, adjacent_SV in zip(self.SVs, self.SVs[1:]):

            ## Compute distance between current and adjacent SV        
            dist = adjacent_SV.beg - SV.end

            ## Add to the list
            distances.append(dist)
        
        return distances

    def shuffle_SVs(self):
        '''
        Randomly shuffle SVs along the chromosome. 
        '''
        ## 1. First randomly shuffle SV list
        random.shuffle(self.SVs)

        ## 2. Reasign SVs at random locations in the chromosome
        # Notes:
        #   - Randomization does not allow overlapping SVs to be consistent with MosaiCatcher (not able to call overlapping SVs)
        ## Initiate variables
        targetIntervals = [(0, self.end)]
        intervalsSV = [] # List containing genomic intervals affected by SVs
        shuffledSVs = []

        ## Reasign one SV at a time
        for SV in self.SVs:

            ## Discard target intervals where the SV will not fit 
            targetIntervals = [interval for interval in targetIntervals if (interval[1] - interval[0]) >= (SV.len() + 1)]

            # Remove SV if does not fit within any interval
            if not targetIntervals:
                continue

            ## Redefine target interval coordinates to make sure the SV will fit within
            targetIntervals = [(interval[0], interval[1] - SV.len()) for interval in targetIntervals]

            ## Map target intervals into 1-0 space
            intervals, mappedIntervals = gRanges.mapRanges(targetIntervals) # Here introduce fix to not get overlapping mapped intervals

            ## Randomly pick one interval
            randomNb = random.uniform(0, 1)
            index = gRanges.overlap_multiranges(randomNb, mappedIntervals)[1][0]
            targetInterval = intervals[index]

            ## Randomly pick one position within the interval
            randomPos = random.randrange(targetInterval[0], targetInterval[1])
            
            ## Relocate SV at new random position
            length = SV.len()
            SV_random = copy.deepcopy(SV)        
            SV_random.beg = randomPos
            SV_random.end = randomPos + length

            ## Add relocated SV to the list containing shuffled SVs
            shuffledSVs.append(SV_random)

            ## Add SV coordinates to the list
            intervalsSV.append((SV_random.beg, SV_random.end))

            ## Redefine target intervals (== genomic regions not affected by SVs)
            targetIntervals = self.complementary_intervals(intervalsSV)

        ## 3. Replace original by shuffled SVs
        self.SVs = shuffledSVs

    def complementary_intervals(self, intervals):
        '''
        Compute set of complementary intervals to input intervals to span the complete chromosome sequence. I.E:
        
        Input intervals:                <—————a—————>		  <—b—>		    			 
        Complementary intervals: *<—-A—->			<————B————>	  <————C———>*
                             chrom_beg                                       chrom_end
        Input:
            1. intervals: list of tuples. Each tuple will correspond to a interval with start and end coordinates. 
                          Note: intervals are asumed not to overlap. If we want to make this more general
                                we should apply an extra step at the beginning to merge overlapping intervals

        Output:
            2. compIntervals: list of complementary interval. Each one represented as a tuple with start and end coordinates
        '''
        ## 1. Sort intervals by start position in increasing coordinates ordering
        intervals.sort(key=lambda range: range[0]) 

        ## 2. Initialize list including first complementary interval (A)
        # Input intervals:                <—————a—————>	  			 
        # Complementary interval:  *<—-A—->			
        #                     chrom_beg (0) 
        #  
        if intervals[0][0] - 0 > 0: # Only create interval if longer than 0bp

            compIntervals = [(0, intervals[0][0])]

        ## 3. Create complementary intervals between adjacent input ranges
        # Input intervals:            <—————a—————>		    <—b—>		    			 
        # Complementary intervals: 	              <————B————>	      
        for interval, adjacent in zip(intervals, intervals[1:]):

            if adjacent[0] - interval[1] > 0: # Only create interval if longer than 0bp

                compIntervals.append((interval[1], adjacent[0]))

        ## 4. Create last complementary interval between last input range and the end of the chromosome
        # Input intervals:         <—b—>		    			 
        # Complementary interval:       <————C———>*
        #                                     chrom_end
        if self.end - intervals[-1][1] > 0: # Only create interval if longer than 0bp

            compIntervals.append((intervals[-1][1], self.end))

        return compIntervals

    def SV_clustering(self, minNbSV, nbRandom):
        '''
        Assess SV clustering. SVs are not randomly distributed along the chromosome. 
        Lower SV distance than expected by chance

        Input:
            1. minNbSV: minimum number of SVs in the chrom to assess clustering
            2. nbRandom: number of randomizations for computing random distribution

        Output: 
            Update three attributes:
            1. clustering: boolean. Cluster with SV clustering (True) or not (False)
            2. c_score: SV clustering score
            3. c_pvalue: p-value resulting from Kolmogórov-Smirnov test
        '''
        ## 1. Abort if chrom does not have enough number of SVs
        if self.nbSV() < minNbSV:
            self.clustering = False 
            self.c_score = None
            self.c_pvalue = None
            return

        ## 2. Compute observed distance distribution
        observed = self.SV_distance()

        ## 3. Compute expected distance distribution 
        randomizations = []

        ## 3.1 Apply multiple randomizations 
        for i in range(nbRandom):

            ## Create copy of the chromosome
            chromosome = copy.deepcopy(self)

            ## Shuffle SVs 
            chromosome.shuffle_SVs()
            
            ## Obtain distance distribution
            randomization = chromosome.SV_distance()

            ## Add random distribution to the list
            randomizations.append(randomization)
        
        ## 3.2 Merge randomizations
        expected = list(itertools.chain.from_iterable(randomizations))

        ## 4. Compare observed with expected distribution (randomized)
        statistic, self.c_pvalue = stats.ks_2samp(observed, expected)

        if self.c_pvalue < 0.05:
            self.clustering = True
        
        else:
            self.clustering = False

        ## 5. Compute clustering score
        pseudocount = 0.0000001
        ratio = (np.median(observed) + pseudocount) / (np.median(expected) + pseudocount)
        self.c_score = math.log(ratio, 2) * -1

    def haplo_bias(self):
        '''
        Assess haplotype bias for the chromosome

        Output: 
            Update three attributes:
            1. haploBias: boolean. Cluster with haplotype bias (True) or not (False) 
            2. h_score: haplotype bias score (fisher exact test odd ratio)
            3. h_pvalue: p-value for haplotype bias          
        '''        
        ## 1. Collect all haplotypes
        haplotypes = []

        # For each SV
        for SV in self.SVs:

            # a) Heterozygous
            if SV.haplo1 in ['h1', 'h2']:
                haplotypes.append(SV.haplo1)
            
            # b) Homozygous (variant in both haplotypes)
            elif SV.haplo1 == 'hom':
                haplotypes = haplotypes + ['h1', 'h2']

            # c) Skip as haplotype info not available

        ## 2. Compute expected distribution (random)
        nbHaplo = len(haplotypes)
        nbH1 = int(nbHaplo / 2)
        nbH2 = int(nbHaplo / 2)
        expected = [nbH1, nbH2]

        ## 3. Compute observed distribution
        haploCounts = collections.Counter({'h1' : 0, 'h2' : 0})
        haploCounts.update(collections.Counter(haplotypes))
        observed = [haploCounts['h1'], haploCounts['h2']]

        ## 4. Compare observed with random through fisher exact test
        self.h_score, self.h_pvalue = stats.fisher_exact([observed, expected])

        if (self.h_pvalue is not None) and (self.h_pvalue < 0.05):
            self.haploBias = True
        
        else:
            self.haploBias = False

    def create_windows_SVs(self, windowSize, minNbSV):
        '''
        Retrieve genomic windows of length == windowSize. Windows will contain 
        the list of SVs falling within.  

        Input:
            1. windowSize: window size
            2. minNbSV: minimum number of SVs in a window to be reported

        Output: 
            1. windows: list of genomic windows containing SVs
        '''
        ## 1. Assign SVs to windows and create them
        windows = {}

        ## Per each SV
        for SV in self.SVs:

            ## Create windows spanned by the SV
            firstWindow = int(SV.beg / windowSize)  
            lastWindow = int(SV.end / windowSize)  

            # For each window
            for windowId in range(firstWindow, lastWindow + 1, 1):
                
                ## Initialize window if needed
                if windowId not in windows:

                    windowBeg = windowId * windowSize
                    windowEnd = (windowId + 1) * windowSize
                    windows[windowId] = genomic_window(self.id, windowBeg, windowEnd)
                    
                ## Add SV to the window
                windows[windowId].add_SV([SV])

        ## 2. Filter out windows without enough number of SVs
        ## 2.1 Create list of windows to remove
        windows2remove = []

        for windowId, window in windows.items():
            if window.nbSV() < minNbSV:
                windows2remove.append(windowId)

        ## 2.2 Remove windows
        for windowId in windows2remove:
            del windows[windowId]
        
        return list(windows.values())

    def group_adjacent(self, maxDist, criteria):
        '''
        Group adjacent SVs

        Input:
            1. maxDist. Maximum distance between two SVs to be considered adjacent 
            2. criteria. List containing criteria used for grouping. Several possible criteria can be used: 
                - 'AF'. Allele frequency
                - 'SVTYPE'. Structural variant type consistency
                - 'HAPLO'. Haplotype consistency

        Output:
            1. groups. Nested list of SV groups. Each inner list will correspond to a SV group. 
        '''

        ## 0. Single SV identified. Return single group
        if (len(self.SVs) == 1):
            return [self.SVs]

        ## 1. Sort SVs along the chromosome
        self.sort_SV()

        ## 2. Group adjacent SVs 
        index2groups = {}
        groupCount = 0

        for index, SV in enumerate(self.SVs[:-1]):
            
            ## 2.1 Pick next SV
            SV_next = self.SVs[index + 1]

            ## 2.2 Compare adjacent SV features
            ## Breakpoint distance  
            dist = SV_next.beg - SV.end

            if (dist <= maxDist):
                DIST_filter = True
            
            else:
                DIST_filter = False

            ## Allele frequency (AF) differences 
            AF_diff = abs(SV.AF - SV_next.AF)
            AF_diff_threshold = 0.30 

            if ('AF' not in criteria) or (AF_diff <= AF_diff_threshold):
                AF_filter = True
            
            else:
                AF_filter = False

            ## SV type consistency
            consistent_SV_type = True if SV.SV_type1 == SV_next.SV_type1 else False

            if ('SVTYPE' not in criteria) or (consistent_SV_type):
                SVTYPE_filter = True
            
            else:
                SVTYPE_filter = False

            ## Haplotype consistency
            consistent_haplotype = True if SV.haplo1 == SV_next.haplo1 else False

            if ('HAPLO' not in criteria) or (consistent_haplotype):
                HAPLO_filter = True
            
            else:
                HAPLO_filter = False

            ## 2.3 Group adjacent SVs passing all the filters
            ## A) Group SVs
            if DIST_filter and AF_filter and SVTYPE_filter and HAPLO_filter:
                
                ## a) First SV already belongs to a group -> Add adjacent SV to the same group
                if index in index2groups:
                    index2groups[index + 1] = groupCount

                ## b) First SV does not belongs to a group -> Create group containing both SVs
                else:
                    groupCount += 1
                    index2groups[index] = groupCount 
                    index2groups[index + 1] = groupCount

            ## B) Do not group SVs 
            else:
                
                ## First SV does not belong to any group -> Create group for it
                if index not in index2groups:
                    groupCount += 1
                    index2groups[index] = groupCount 

                ## Create group containing second SV
                groupCount += 1
                index2groups[index + 1] = groupCount 

        ## Rearrange dict to have groups as keys with the list of indexes per group as value 
        groupIndexes = structures.groupKeysByValue(index2groups)

        ## Collect events composing each group based on the indexes
        groups = []

        # For each group 
        for indexes in groupIndexes.values():

            ## Select SVs composing the group based on the indexes
            SVs_in_group = operator.itemgetter(*indexes)(self.SVs)

            ## Group composed by a single SV -> create single element list
            if not isinstance(SVs_in_group, collections.Iterable):
                SVs_in_group = [SVs_in_group]

            ## Group composed by multiple SVs -> convert tuple to list
            else:
                SVs_in_group = list(SVs_in_group)

            ## Add group to the output list
            groups.append(SVs_in_group)
        
        return groups

    def cluster_SV(self):
        '''
        Group adjacent structural variants (SV) into clusters

        Output: Add clusters to 'SV_clusters' attribute
        '''        
        ## 1. Group adjacent SVs
        criteria = []
        groups = self.group_adjacent(400000, criteria)

        ## 2. Create clusters
        for group in groups:
            cluster = variants.SV_cluster(group)
            self.SV_clusters.append(cluster)

    def merge_SV(self):
        '''
        Fix SV oversegmentation by merging adjacent SVs with enough evidences 
        pointing to they correpond to the same SV (same SV type, haplotype, similar allele freq)

        Output: Updated list of SVs along the chromosomes after merging
        '''
        ## 1. Group adjacent SVs prior merging
        criteria = ['SVTYPE']
        groups = self.group_adjacent(0, criteria)

        ## 2. Do merging of adjacent SV groups  
        mergedSVs = []

        for group in groups:

            ## Merge SVs in group
            mergedSV = variants.merge_SV(group)
            mergedSVs.append(mergedSV)

        ## 3. Replace original SVs by merged SVs
        self.SVs = mergedSVs

        ## 4. Sort merged SVs
        self.sort_SV()

    def telomeric_del(self, arm, maxDist2end):
        '''
        Check if the chromosome has a telomeric deletion at a given end

        Input: 
            1. arm: chromosomal arm to search for BFB ('p' or 'q')
            2. maxDist2end: maximum distance of the deletion to the telomere end

        Output:
            1. telDel: SV object for telomeric deletion. None if not found
        '''
        ## 0. Initialization
        telDel = None

        if not self.SVs:
            return telDel

        ## 1. Sort SVs in increasing coordinates ordering
        self.sort_SV()

        ## 2. Compute distance between telomere and adjacent SV
        # a) p-arm
        if (arm == 'p'):
            SV = self.SVs[0]
            dist2end = SV.beg - self.beg    
        
        # b) q-arm
        elif (arm == 'q'):
            SV = self.SVs[-1]
            dist2end = self.end - SV.end 

        # c) Not supported arm
        else:
            print('[ERROR] Not supported arm')
            return telDel

        ## 3. Assess if telomeric deletion:
        # a) Adjacent SV is a deletion AND
        # b) Deletion within maxDist2end with respect to the telomere
        if (SV.SV_type1 == 'del') and (dist2end <= maxDist2end):
            telDel = SV
            telDel.dist2end = dist2end

        return telDel

    def adjacent_dup(self, arm, maxDist):
        '''
        Search for duplication cluster adjacent to a telomeric deletion

        Input: 
            1. arm: chromosomal arm affected by the deletion ('p' or 'q')
            2. maxDist: maximum distance between adjacent events to group them together

        Output:
            1. cluster: list of dup, idup and complex SVs adjacent to the telomeric deletion. Empty list if not found
        '''
        ## 1. Order SVs relative the affected arm
        SVs = self.SVs if arm == 'p' else self.SVs[::-1]
 
        ## 2. Search for a candidate duplication cluster adjacent to the telomeric deletion
        candidate = []

        for SV, SV_next in zip(SVs, SVs[1:]):

            ## Compute distance between consecutive SVs:
            dist = SV_next.beg - SV.end if arm == 'p' else SV.beg - SV_next.end

            if (dist <= maxDist) and (SV_next.SV_type1 in ['complex', 'dup', 'idup']):
                candidate.append(SV_next)
            
            else:
                break


        ## 3. Filter out candidate cluster if it does not contain any dup or idup event
        SV_types = list(itertools.chain.from_iterable([[SV.SV_type1, SV.SV_type2] for SV in candidate])) 

        if ('dup' in SV_types) or ('idup' in SV_types):
            cluster = candidate
        else:
            cluster = []
        
        return cluster
                  
    def search4BFB(self, arm):
        '''
        Search for Breakage Fusion Bridge (BFB) cycle events 

        Input:
            1. arm: chromosomal arm to search for BFB ('p' or 'q')

        Output: 
            1. BFB: BFB object. None if no BFB found
        '''
        ## 1. Assess if chromosome has a telomeric deletion        
        telDel = self.telomeric_del(arm, 500000)

        ## 2. Search for adjacent duplications
        dups = self.adjacent_dup(arm, 7000000)

        ## 3. Create BFB if both telomeric deletion and adjacent duplication found
        if (telDel is not None) and (dups):
            BFB = variants.BFB(self.id, arm, telDel, dups, self.cellId)
            
        else:
            BFB = None
        
        return BFB


    def is_chromothriptic(self):
        '''
        Assess if the chromomosome has undergone chromothripsis

        Output: 
            1. chromothriptic: boolean. Chromothriptic chromosome (True) or not (False) 
        '''
        ## Make chromothripsis call if:
        # 1. Chromosome has enrichment on number of SVs AND 
        # 2. Chromosome has SV clustering (not random distribution) AND
        # 3. Chromosome has bias toward one haplotype
        if self.enrichment and self.clustering and self.haploBias:
            chromothriptic = True
        else:
            chromothriptic = False
        
        return chromothriptic

class genomic_window():
    '''
    Class representing a genomic windows
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end):
        '''
        Initialize class instance
        '''
        genomic_window.number += 1 # Update instances counter
        self.id = 'WINDOW_' + str(genomic_window.number)
        self.ref = ref
        self.beg = beg
        self.end = end
        self.SVs = []

    def add_SV(self, SVs):
        '''
        Add a set of SV_call events to the genomic window

        Input:
            1. SVs: list of SV_call objects

        Output: 
            Update SVs attribute with new SV_call events
        '''
        self.SVs = self.SVs + SVs

    def sort_SV(self):
        '''
        Sort SVs along the chromosome by their starting position (increasing)
        '''            
        self.SVs.sort(key=lambda SV: SV.beg, reverse=False)

    def nbSV(self):
        '''
        Compute the number of SVs in the chromosome

        Output:
            1. nbSVs: integer with the number of SVs
        '''
        return len(self.SVs)

    def median_SV_distance(self):
        '''
        Compute the median distance between adjacent SVs along the interval

        Output:
            1. medianDist: median distance between adjacent SVs
        '''
        ## 1. Retrieve distance distribution
        distances = self.SV_distance()

        ## 2. Compute median distance
        medianDist = statistics.median(distances)

        return medianDist

    def SV_distance(self):
        '''
        Compute the distances between adjacent SVs along the chromosome

        Output:
            1. distances: list containing the distance between the begin and end position of adjacent SVs
        '''
        ## 1. Make sure that SVs are sorted by their start position
        self.sort_SV()

        ## 2. Compute distances between adjacent SVs
        distances = []

        ## For each SV and adjacent SV pair
        for SV, adjacent_SV in zip(self.SVs, self.SVs[1:]):

            ## Compute distance between current and adjacent SV        
            dist = adjacent_SV.beg - SV.end

            ## Add to the list
            distances.append(dist)
        
        return distances     
