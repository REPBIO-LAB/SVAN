'''
Module 'gRanges' - Contains functions to do operations with genomic ranger or coordinates (i.e. overlaps...)
'''

## DEPENDENCIES ##

## FUNCTIONS ##
def makeGenomicBins(refLengths, binSize, targetRefs):
    '''
    '''
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


def rangeList2dict(ranges):
    '''
    Organize list of ranges into a dictionary  

    Input:
        1. ranges: list of ranges. Each list item corresponds to a tuple (ref, beg, end)

    Output:
        1. rangesDict: Dictionary with reference ids as keys and the list of ranges on each reference as values
    '''    
    rangesDict = {}

    # For each interval
    for interval in ranges:
        ref, beg, end = interval

        # Initialize interval list at reference
        if ref not in rangesDict:
            rangesDict[ref] = []

        # Add interval to the list
        rangesDict[ref].append([beg, end])
    
    return rangesDict

def overlap(begA, endA, begB, endB):
    '''
    Check if two ranges overlap    

    Input:
        1. begA: interval A begin position
        2. endA: interval A end position
        3. begB: interval B begin position
        4. endB: interval B end position
    Output:
        1. boolean: True (overlap) and False (no overlap)
        2. overlapLen: number of overlapping base pairs
    '''    
    maxBeg = max([begA, begB])
    minEnd = min([endA, endB])

    # a) Overlap
    if maxBeg <= minEnd:
        boolean = True
        overlapLen = minEnd - maxBeg + 1

    # d) No overlap
    else:
        boolean = False
        overlapLen = 0

    return boolean, overlapLen

def overlap_extended(begA, endA, begB, endB):
    '''
    Check if two ranges overlap. Report also the overlapping coordinates for each range

    Input:
        1. begA: interval A begin position
        2. endA: interval A end position
        3. begB: interval B begin position
        4. endB: interval B end position
    Output:
        1. boolean: True (overlap) and False (no overlap)
        2. overlapLen: number of overlapping base pairs
        3. coord: tuple with overlapping coordinates (None if no overlap found)
    '''    
    ## a) No overlap
    # <-----A----->endA   begB<-----B-----> OR
    # <-----B----->endB   begA<-----A----->
    if (begB > endA) or (begA > endB):
        boolean = False
        overlapLen = 0
        coord = None

    ## b) Range A within B
    #     begA<-----A----->endA
    # begB<---------B--------->endB
    elif (begA >= begB) and (endA <= endB):
        boolean = True
        overlapLen = endA - begA 
        coord = (begA, endA)

    ## c) Range B within A
    # begA<---------A--------->endA
    #     begB<-----B----->endB
    elif (begB >= begA) and (endB <= endA):
        boolean = True
        overlapLen = endB - begB 
        coord = (begB, endB)

    ## d) Partial overlap (A first)
    #   begA<-----A----->endA   
    #           begB<-----B----->   
    elif (begB >= begA) and (begB <= endA):
        boolean = True
        overlapLen = endA - begB 
        coord = (begB, endA)

    ## e) Partial overlap (B first)
    #          begA<-----A----->
    #   begB<-----B----->endB
    elif (begA >= begB) and (begA <= endB):
        boolean = True
        overlapLen = endB - begA 
        coord = (begA, endB)

    return boolean, overlapLen, coord

def rcplOverlap(begA, endA, begB, endB, minPercOverlap):
    '''
    Check if two ranges overlap with a minimum percentage of reciprocal overlap

    Input:
        1. begA: interval A begin position
        2. endA: interval A end position
        3. begB: interval B begin position
        4. endB: interval B end position
        5. minPercOverlap: minimum percentage of reciprocal overlap between A and B

    Output:
        1. boolean: True (overlap) and False (no overlap)
        2. overlapLen: number of overlapping base pairs
    '''    
    overlapLen = overlap(begA, endA, begB, endB)[1]
    lenA = (endA - begA)
    lenB = (endB - begB)
    percA = (overlapLen / lenA) * 100 if lenA > 0 else 0
    percB = (overlapLen / lenB) * 100 if lenB > 0 else 0

    # a) Both ranges overlap with a minimum X percentage of reciprocal overlap
    if (percA >= minPercOverlap) and (percB >= minPercOverlap):
        boolean = True

    # b) No overlap
    else:
        boolean = False

    return boolean, overlapLen


def complementary(begA, endA, begB, endB, maxDist, maxPercOverlap):
    '''
    Check if two intervals are complementary. Two ranges are considered complementary when:
    1) Ranges located less than X bp of distance AND
    2) Overlapping region does not span X% or more of any of the intervals

    Input:
        1. begA: begin position of range A
        2. endA: end position of range A
        3. begB: begin position of range B
        4. endB: end position of range B
        5. maxDist: maximum distance between both ranges 
        6. maxPercOverlap: maximum percentage of overlap between ranges

    Output:
        1. boolean: True (intervals are complementary) or False (intervals not complementary)
        2. orientation: 'LEFT' (B on the left of A), 'RIGHT' (B on the right of A) or None (not complementary)
    '''  
    ## 1. Redefine intervals by adding the maximum distance to their begin and end coordinates
    newBegA = int(begA) - maxDist
    newEndA = int(endA) + maxDist
    newBegB = int(begB) - maxDist
    newEndB = int(endB) + maxDist

    ## 2. Assess if redefined intervals do overlap
    boolean = overlap(newBegA, newEndA, newBegB, newEndB)[0]
    
    # A) No overlap
    if not boolean:
        orientation = None
    
    # B) Overlap 
    else:

        # 3. Discard those cases with an overlapping region spanning X% or more of at least one of the intervals
        # --------A---------
        #              <---> Overlapping region
        #              -------B-------
        overlapLen = overlap(begA, endA, begB, endB)[1]   # Compute real degree of overlap

        # For each interval, compute its % overlapping the other interval
        lenA = (endA - begA) + 1
        lenB = (endB - begB) + 1
        percA = (overlapLen / lenA) * 100
        percB = (overlapLen / lenB) * 100

        # a) One of the intervals overlaps the other by more than the maximum allowed % cutoff 
        if (percA > maxPercOverlap) or (percB > maxPercOverlap):
            boolean = False
            orientation = None

        # b) Overlapping region does not span X% or more of any of the intervals
        else:
            boolean = True

            ## Determine complementariety orientation (use A as reference)
            # a) B located on the left of A
            #                   begA <------A------> 
            #  begB <------B------>
            if begA > begB:
                orientation = 'LEFT'

            # b) B located on the right of A
            #  begA <------A------>
            #                  begB <------B------> 
            else:
                orientation = 'RIGHT'

    return boolean, orientation