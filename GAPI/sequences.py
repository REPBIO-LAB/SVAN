'''
Module 'sequences' - Contains functions for the manipulation and extracting information from nucleotidic sequences
'''

## DEPENDENCIES ##
# External
import os
import subprocess

# Internal
from GAPI import formats
from GAPI import log
from GAPI import unix

## FUNCTIONS ##
def rev_complement(seq):
    '''
    Make the reverse complementary of a dna sequence

    Input:
        1. seq: DNA sequence
    Output:
        1. revComplementSeq: Reverse complementary of input DNA sequence
    '''
    baseComplementDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', ',': 'N'}
    seq = seq.upper()
    revSeq = seq[::-1] # Make reverse sequence
    letters = list(revSeq)
    letters = [baseComplementDict[base] for base in letters]
    revComplementSeq = ''.join(letters) # Make complement of reverse sequence

    return revComplementSeq


def baseComposition(seq):
    '''
    Compute mononucleotide frequencies for an input sequence

    Input:
        1. seq: DNA sequence
    Output:
        1. baseCounts: Dictionary containing counts per nucleotidic base for the input sequence 
        1. basePercs: Dictionary containing percentajes per nucleotidic base for the input sequence
    '''
    ## 1. Compute counts
    baseCounts = {}
    baseCounts['A'] = seq.count('A') 
    baseCounts['G'] = seq.count('G') 
    baseCounts['C'] = seq.count('C') 
    baseCounts['T'] = seq.count('T') 
    baseCounts['N'] = seq.count('N')
    baseCounts['total'] = baseCounts['A'] + baseCounts['G'] + baseCounts['C'] + baseCounts['T'] + baseCounts['N']

    ## 2. Compute percentages
    basePercs = {}
    seqLen = len(seq)

    for base in baseCounts:
        perc = baseCounts[base] / seqLen * 100
        basePercs[base] = perc

    return baseCounts, basePercs


def find_monomers(seq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity):
    '''
    Find nucleotide monomers in an input sequence

    Input:
        1. seq: DNA sequence
        2. targetMonomer: Type of monomer to search for
        3. windowSize: Window size used to scan the sequence
        4. maxWindowDist: Maximum number of windows without the monomer between two monomeric windows to group them in the same cluster
        5. minMonomerSize: Minimum size for calling a monomer
        6. minPurity: minimum % of bases corresponding to the target monomer in a window or monomer cluster

    Output:
        1. monomers: List containing the monomer instances identified 
    '''
    seqLen = len(seq)

    windowDist = 0
    monomers = []
    noMonomeric = []
    indexes = []

    ## 1. Split seq into consecutive slices 
    # ---ws--->---ws--->---ws--->---ws--->
    slices = [(pos, pos + windowSize, seq[pos:pos + windowSize]) for pos in range(0, seqLen, windowSize)]

    ## 2. Search for candidate monomers by iterating over the slices
    for index, value in enumerate(slices):
        
        indexes.append(index)

        # Collect slide info
        beg, end, sliceSeq = value
        sliceLen = len(sliceSeq)

        # Compute % of target monomer in the slice 
        baseCounts, basePercs = baseComposition(sliceSeq)
        percMonomer = float(baseCounts[targetMonomer])/sliceLen*100
        
        ## A) Slice correspond to a monomer (---ws--- == AAA.../TTT.../...) if:
        # 1. sliceLen >= windowSize AND 
        # 2. purity >= minPurity
        if (sliceLen >= windowSize) and (percMonomer >= minPurity):

            ## a) CREATE new monomer object if:
            #  a.1) No monomer already identified OR
            #  a.2) Distance from new monomer to the previously reported monomer > maxWindowDist 
            if (not monomers) or (windowDist > maxWindowDist):

                ## Create monomer
                monomerObj = monomer(beg, sliceSeq, index)
                monomers.append(monomerObj)

            ## b) EXTENSION if monomer within maxWindowDist to previous identified monomer 
            else:
                previousMonomer = monomers[-1]
                
                ## b.1 Extend previously identified monomer with all the not monomeric slices till the new monomer 
                for noMonomericSlice in noMonomeric:
                    noMonomericSliceSeq, noMonomericIndex = noMonomericSlice[2:]
                    previousMonomer.add(noMonomericSliceSeq, noMonomericIndex)

                ## b.2 Extend previously identified monomer with the new monomer         
                previousMonomer.add(sliceSeq, index)

            ## Initiate distance counter and no monomeric list
            windowDist = 0
            noMonomeric = []

        ## B) Slide do not correspond to a monomer
        else:
            windowDist += 1
            noMonomeric.append([beg, end, sliceSeq, index])
    
    ## 3. Refine monomers by attempting to extend both ends 
    # [ACGAA] AAAAAAAAAAAAAAAAAAAA [AAAGG]
    # Look to the slices flanking the monomer and extend till no A match is found
    # In the example the monomer would be extended as follows:
    # [ACG|AA]AAAAAAAAAAAAAAAAAAAA[AAA|GG]
    # '|' indicates the end of the extension 
    
    # For each monomer
    for monomerObj in monomers:
    
        ## 3.1 Left extension
        # [ACG|AA]AAAAAAAAAAAAAAAAAAAA
        # '|' indicates the end of the extension 
        index = monomerObj.indexes[0] - 1 

        # flanking slice available
        if index in indexes:

            # Pick sequence
            sliceSeq = slices[index][2]

            # Parse sequence backward one nucleotide at a time
            # [ACG|AA] *Monomer* 
            # <-<-<-<-
            for nucleotide in reversed(sliceSeq):

                # a) Extend monomer if nucleotide match
                if nucleotide == targetMonomer:
                    
                    monomerObj.seq = nucleotide + monomerObj.seq
                    monomerObj.beg -= 1

                    # Add slice index to the list
                    if index not in monomerObj.indexes:
                        monomerObj.indexes.insert(0, index)
 
                
                # b) Stop extension once no nucleotide match is found
                else:
                    break

        ## 3.2 Right extension
        # AAAAAAAAAAAAAAAAAAAA[AAA|GG]
        # '|' indicates the end of the extension 
        index = monomerObj.indexes[-1] + 1 
        
        # flanking slice available
        if index in indexes:

            # Pick sequence
            sliceSeq = slices[index][2]

            # Parse sequence forward one nucleotide at a time
            # *Monomer* [AAA|GG]
            #           ->->->->
            for nucleotide in sliceSeq:

                # a) Extend monomer if nucleotide match
                if nucleotide == targetMonomer:
                    
                    monomerObj.seq = monomerObj.seq + nucleotide 
                    monomerObj.end += 1

                    # Add slice index to the list
                    if index not in monomerObj.indexes:
                        monomerObj.indexes.append(index) 

                # b) Stop extension once no nucleotide match is found
                else:
                    break            

    ## 4. Filter candidate monomers
    filteredMonomers = []

    # For each monomer
    for monomerObj in monomers:
        monomerSize = len(monomerObj.seq)
        
        ## Select those monomers with:
        # 1. size >= minMonomerSize 
        if (monomerSize >= minMonomerSize):
            filteredMonomers.append(monomerObj)
        
    return monomers

def find_monomers_backward(seq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity):
    '''
    Find nucleotide monomers in an input sequence

    Input:
        1. seq: DNA sequence
        2. targetMonomer: Type of monomer to search for
        3. windowSize: Window size used to scan the sequence
        4. maxWindowDist: Maximum number of windows without the monomer between two monomeric windows to group them in the same cluster
        5. minMonomerSize: Minimum size for calling a monomer
        6. minPurity: minimum % of bases corresponding to the target monomer in a window or monomer cluster

    Output:
        1. monomers: List containing the monomer instances identified 
    '''
    seqLen = len(seq)

    windowDist = 0
    monomers = []
    noMonomeric = []
    indexes = []

    ## 1. Split seq into consecutive slices 
    # <---ws---<----ws---<---ws---<---ws---
    slices_rev = [(pos, pos - windowSize, seq[pos - windowSize:pos][::-1]) if ((pos - windowSize) >= 0) else (0, pos, seq[0:pos]) for pos in range(seqLen, 0, (windowSize * -1))]

    ## 2. Search for candidate monomers by iterating over the slices
    for index, value in enumerate(slices_rev):
    
        indexes.append(index)

        # Collect slide info
        end, beg, sliceSeq = value
        sliceLen = len(sliceSeq)

        # Compute % of target monomer in the slice 
        baseCounts, basePercs = baseComposition(sliceSeq)
        percMonomer = float(baseCounts[targetMonomer])/sliceLen*100
        
        ## A) Slice correspond to a monomer (---ws--- == AAA.../TTT.../...) if:
        # 1. sliceLen >= windowSize AND 
        # 2. purity >= minPurity
        if (sliceLen >= windowSize) and (percMonomer >= minPurity):

            ## a) CREATE new monomer object if:
            #  a.1) No monomer already identified OR
            #  a.2) Distance from new monomer to the previously reported monomer > maxWindowDist 
            if (not monomers) or (windowDist > maxWindowDist):

                ## Create monomer
                monomerObj = monomer(beg, sliceSeq, index)
                monomers.append(monomerObj)

            ## b) EXTENSION if monomer within maxWindowDist to previous identified monomer 
            else:
                previousMonomer = monomers[-1]
                
                ## b.1 Extend previously identified monomer with all the not monomeric slices till the new monomer 
                for noMonomericSlice in noMonomeric:
                    noMonomericSliceSeq, noMonomericIndex = noMonomericSlice[2:]
                    previousMonomer.add_backward(noMonomericSliceSeq, noMonomericIndex)

                ## b.2 Extend previously identified monomer with the new monomer         
                previousMonomer.add_backward(sliceSeq, index)

            ## Initiate distance counter and no monomeric list
            windowDist = 0
            noMonomeric = []

        ## B) Slide do not correspond to a monomer
        else:
            windowDist += 1
            noMonomeric.append([beg, end, sliceSeq, index])
    
    ## 3. Refine monomers by attempting to extend both ends 
    # [ACGAA] AAAAAAAAAAAAAAAAAAAA [AAAGG]
    # Look to the slices flanking the monomer and extend till no A match is found
    # In the example the monomer would be extended as follows:
    # [ACG|AA]AAAAAAAAAAAAAAAAAAAA[AAA|GG]
    # '|' indicates the end of the extension 
    
    # For each monomer
    for monomerObj in monomers:

        ## 3.1 Left extension
        # [ACG|AA]AAAAAAAAAAAAAAAAAAAA
        # '|' indicates the end of the extension 
        index = monomerObj.indexes[0] - 1 

        # flanking slice available
        if index in indexes:

            # Pick sequence
            sliceSeq = slices_rev[index][2]

            # Parse sequence backward one nucleotide at a time
            # [ACG|AA] *Monomer* 
            # <-<-<-<-
            for nucleotide in reversed(sliceSeq):

                # a) Extend monomer if nucleotide match
                if nucleotide == targetMonomer:
                    
                    monomerObj.seq = nucleotide + monomerObj.seq
                    monomerObj.end += 1

                    # Add slice index to the list
                    if index not in monomerObj.indexes:
                        monomerObj.indexes.insert(0, index)
 
                # b) Stop extension once no nucleotide match is found
                else:
                    break

        ## 3.2 Right extension
        # AAAAAAAAAAAAAAAAAAAA[AAA|GG]
        # '|' indicates the end of the extension 
        index = monomerObj.indexes[-1] + 1 
        
        # flanking slice available
        if index in indexes:

            # Pick sequence
            sliceSeq = slices_rev[index][2]

            # Parse sequence forward one nucleotide at a time
            # *Monomer* [AAA|GG]
            #           ->->->->
            for nucleotide in sliceSeq:

                # a) Extend monomer if nucleotide match
                if nucleotide == targetMonomer:
                    
                    monomerObj.seq = monomerObj.seq + nucleotide 
                    monomerObj.beg -= 1

                    # Add slice index to the list
                    if index not in monomerObj.indexes:
                        monomerObj.indexes.append(index) 

                # b) Stop extension once no nucleotide match is found
                else:
                    break    

    ## 4. Filter candidate monomers 
    filteredMonomers = []

    # For each monomer
    for monomerObj in monomers:
        monomerSize = len(monomerObj.seq)
        
        ## Select those monomers with size >= minMonomerSize 
        if (monomerSize >= minMonomerSize):

            ## Set right sequence orientationorientation
            monomerObj.seq = monomerObj.seq[::-1]

            ## Add to the filtered list
            filteredMonomers.append(monomerObj)
    
    return filteredMonomers[::-1]

def filter_internal_monomers(monomers, targetSeq, maxDist2Ends, minInternalMonomerSize):
    '''
    Filter monomers by requesting a larger size for internal monomers

    Input:
        1. monomers: 
        2. targetSeq:
        3. maxDist2Ends: 
        4. minInternalMonomerSize: 

    Output:
        1. filteredMonomers: 
    '''
    ## An internal monomer is located at more than Xbp from the inserted sequence ends
    filteredMonomers = []

    ## For each input monomer
    for monomer in monomers:
        dist2Beg = monomer.beg
        dist2End = len(targetSeq) - monomer.end

        baseCounts, basePercs = baseComposition(monomer.seq)
        percMonomer = float(baseCounts['A'])/monomer.length()*100

        ## a) Select external monomers located at less than or equal to X bp from insert ends 
        if (dist2Beg <= maxDist2Ends) or (dist2End <= maxDist2Ends):
            filteredMonomers.append(monomer)
                            
        ## b) Select internal monomers longer or equal than X bp
        elif (monomer.length() > minInternalMonomerSize) and (percMonomer > 90):
            filteredMonomers.append(monomer)
            
        ## c) Discard internal monomers smaller than X bp

    return filteredMonomers

def filter_no3prime_monomers(monomers, targetSeq, maxDist2end, minInternalMonomerSize):
    '''
    Filter monomers by requesting a larger size for monomers not at the 3' end

    Input:
        1. monomers: 
        2. targetSeq:
        3. maxDist2end: 
        4. minInternalMonomerSize: 

    Output:
        1. filteredMonomers: 
    '''
    minDist2beg = maxDist2end

    ## A no 3' monomer is located at more than Xbp from the inserted sequence 3' end
    filteredMonomers = []

    ## For each input monomer
    for monomer in monomers:
        dist2End = len(targetSeq) - monomer.end
        dist2beg = monomer.beg

        baseCounts, basePercs = baseComposition(monomer.seq)
        percMonomer = float(baseCounts['A'])/monomer.length()*100

        ## a) Select external monomers located at less than or equal to X bp from insert ends 
        if (dist2End <= maxDist2end):
            filteredMonomers.append(monomer)
                            
        ## b) Select internal monomers longer or equal than X bp
        elif (dist2beg >= minDist2beg) and (monomer.length() > minInternalMonomerSize) and (percMonomer > 90):
            filteredMonomers.append(monomer)
            
        ## c) Discard internal monomers smaller than X bp

    return filteredMonomers

## CLASSES ##
class monomer():
    '''
    '''
    def __init__(self, beg, seq, index):
        '''
        Initialize empty class instance
        '''
        self.beg = beg
        self.end = self.beg + len(seq)
        self.seq = seq
        self.indexes = [index]

    def add(self, seq, index):
        '''
        Extend s by adding a new piece of sequence at its end

        Input:
            1. seq: piece of sequence to add at monomer end
            2. index: slice index
        '''
        self.end = self.end + len(seq)
        self.seq = self.seq + seq
        self.indexes.append(index)

    def add_backward(self, seq, index):
        '''
        Extend s by adding a new piece of sequence at its end

        Input:
            1. seq: piece of sequence to add at monomer end
            2. index: slice index
        '''
        self.beg = self.beg - len(seq)
        self.seq = self.seq + seq
        self.indexes = self.indexes + [index]

    def length(self):
        '''
        Return monomer length
        '''
        return len(self.seq)

def create_targeted_fasta(targetIntervalList, reference, outDir):
    '''
    Extract regions of interest from a fasta file.
    
    Input:
        1. targetIntervalList: Reference genome list of intervals to be extracted. The intervals must be provided as chr:beg-end.
        2. reference: Path to fasta file. An index of the reference generated with samtools faidx must be located in the same directory
        3. outDir: Output directory

    Output:
        1. target: Path to fasta file with sequences extarcted from intervals.
    '''
    ## 0. Create logs directory
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Extract the reference target regions 
    target = outDir + '/targetRegions.fa'
    err = open(logDir + '/target.err', 'w')
    targetRegionsPath = outDir + '/targetRegions.txt'
    targetRegions = open(targetRegionsPath, 'w')
    
    for targetInterval in targetIntervalList:
        targetRegions.write(targetInterval + '\n')
    targetRegions.close()

    command = 'samtools faidx ' + reference + ' -r ' + targetRegionsPath + ' -o ' + target
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'TARGET'
        msg = 'Extraction of reference target region failed' 
        log.step(step, msg)
        return None

    # TODO: remove targetRegionsPath file

    return target

def komplexityFilter(komplexityThreshold, inFasta, outFasta, outDir):
    '''
    Filter fasta file using komplexity tool
    Input:
        1. komplexityThreshold: Complexity threshold filter.
        2. inFasta: input FASTA file name
        3. outFasta: output FASTA file name
        4. outDir: input AND output directory (it must be the same)
    Output:
        1. allFastas: Filteres FASTA file complete path.
    '''

    # Set input an output files
    allFastas_all = outDir + '/' + inFasta
    allFastas = outDir + '/' + outFasta

    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    command = 'kz --filter --threshold ' + str(komplexityThreshold) + ' --fasta < ' + allFastas_all + ' > ' + allFastas
    err = open(logDir + '/komplexity.err', 'w') 
    status = subprocess.call(command, stderr=err, shell=True)
    if status != 0:
        step = 'KOMPLEXITY'
        msg = 'Komplexity filter failed. PID: ' + str(os.getpid())
        log.step(step, msg)
    
    return allFastas