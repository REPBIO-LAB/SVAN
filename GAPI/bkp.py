'''
Module to solve the bkp 
'''

## External
import os
import subprocess
import operator
import re
import pysam
import itertools

## Internal
from GAPI import databases
from GAPI import formats
from GAPI import log
from GAPI import unix
from GAPI import sequences
from GAPI import assembly
from GAPI import alignment
from GAPI import events
from GAPI import bamtools
from GAPI import clusters

def bkp(metaclusters):
    '''
    Determine breakpoints
    
    Input: 
    1. List of metaclusters
    
    Output:
    There is no output, but refLeftBkp and refRightBkp attributes are filled
    '''
    for metacluster in metaclusters:
        
        leftBkps = []
        rightBkps = []

        # Choose bkp with highest number of clipping events
        for event in metacluster.events:
            
            if event.type == 'CLIPPING':
                
                if event.clippedSide == 'left':
                    leftBkps.append(event.beg)
                    
                elif event.clippedSide == 'right':
                    rightBkps.append(event.beg)

        # Fill refLeftBkp and refRightBkp attributes
        if len(leftBkps) > 0:
            leftBkp = max(set(leftBkps), key=leftBkps.count)
            metacluster.refLeftBkp = leftBkp
            
        if len(rightBkps) > 0:
            rightBkp = max(set(rightBkps), key=rightBkps.count)
            metacluster.refRightBkp = rightBkp


def bkp_fromClusters(metaclusters):
    '''
    Determine breakpoints
    
    Input: 
    1. List of metaclusters
    
    Output:
    There is no output, but refLeftBkp and refRightBkp attributes are filled
    '''
    for metacluster in metaclusters:
        
        leftBkps = []
        rightBkps = []
        
        ## A) Create subclusters
        subclusters = metacluster.create_subclusters()
        
        ## B) Collect read clipping boundaries
        # if there is a LEFT-CLIPPING cluster:
        if 'LEFT-CLIPPING' in subclusters.keys():
            leftBkps = [event.beg for event in subclusters['LEFT-CLIPPING'].events]
        
        # elif there is a LEFT-CLIPPING cluster
        elif 'RIGHT-CLIPPING' in subclusters.keys():
            rightBkps = [event.beg for event in subclusters['RIGHT-CLIPPING'].events]
        
        ## C) Set bkps
        # if right and left bkps found
        if leftBkps and rightBkps:
            
            leftBkp = max(set(leftBkps), key=leftBkps.count)
            rightBkp = max(set(rightBkps), key=rightBkps.count)
            
            bkp = max(set(rightBkps+leftBkps), key=rightBkps.count)
            
            # if there is less than 20 bp between bkps
            if abs(rightBkp - leftBkp) < 20:
                metacluster.refLeftBkp, metacluster.refRightBkp = leftBkp, rightBkp if leftBkp < rightBkp else bkp, bkp
        
        # if just left bkp found
        elif leftBkps:
            leftBkp = max(set(leftBkps), key=leftBkps.count)
            metacluster.refLeftBkp = leftBkp
        
        # if just right bkp found
        elif rightBkps:
            rightBkp = max(set(rightBkps), key=rightBkps.count)
            metacluster.refRightBkp = rightBkp

                
def analyzeMetaclusters(metaclusters, confDict, bam, normalBam, mode, outDir, binId, identDbPath):
    '''
    Four main steps
        a. Add supporting clipping events to discordant metacluster.
        b. Determine metacluster bkps.
        c. Reconstruct bkp sequence.
        d. Determine inserted sequence bkp.
    
    Input
        1. metaclusters: List of metaclusters
        2. confDict
        3. bam
        4. normalBam
        5. mode
        6. outDir
        7. binId
        8. identDbPath
    
    Output
        Doesn't return anything, just fill some metaclusters attributes:
        1. metacluster.refRightBkp
        2. metacluster.refLeftBkp
        3. metacluster.rightSeq
        4. metacluster.leftSeq
        5. metacluster.intRightBkp
        6. metacluster.intLeftBkp
    '''

    # Initiliaze dictionary: metaclustersWODiscClip[metacluster] = clippings -> key: metacluster without discordant clipping events; value: list of candidate clippings events to add.
    metaclustersWODiscClip = {}

    # Make bkp directory
    bkpDir = outDir + '/BKP'
    unix.mkdir(bkpDir)

    # a. Add supporting clipping to discordant metacluster.
    
    for metacluster in metaclusters:

        # Make bkp directory
        specBkpDir = bkpDir + '/' + metacluster.ref + '_' + str(metacluster.beg) + '_' + str(metacluster.end)
        unix.mkdir(specBkpDir)

        if metacluster.orientation != 'RECIPROCAL':
            clippingEventsToAdd, discClip = supportingCLIPPING(metacluster, 150, confDict, bam, normalBam, mode, metacluster.orientation)
            # TODO SR: For the moment, no BLAT id performed for MEs insertions.
            # Add clippings if there are discodant clippings.
            if discClip == True or not 'VIRUS' in confDict['targetINT2Search'] or not metacluster.identity:
                metacluster.addEvents(clippingEventsToAdd)
                # Set bkp clipping type
                if discClip == True:
                    if metacluster.orientation == 'PLUS':
                        metacluster.rightClipType = 'DISC'
                    elif metacluster.orientation == 'MINUS':
                        metacluster.leftClipType = 'DISC'

                else:
                    if metacluster.orientation == 'PLUS':
                        metacluster.rightClipType = 'REG'
                    elif metacluster.orientation == 'MINUS':
                        metacluster.leftClipType = 'REG'

            # If there are not discordant clippings, but metacluster has identity, return metacluster and clippings list in order to perform BLAT clippings search
            elif clippingEventsToAdd and not discClip and metacluster.identity:
                metaclustersWODiscClip[metacluster] = clippingEventsToAdd
            # If there are not discordant clippings and the metacluster has no identity, dont add clippings.

        elif metacluster.orientation == 'RECIPROCAL':
            clippingRightEventsToAdd, discClip = supportingCLIPPING(metacluster, 150, confDict, bam, normalBam, mode, 'PLUS')
            # Add clippings if there are discodant clippings.
            if discClip == True or not 'VIRUS' in confDict['targetINT2Search'] or not metacluster.identity:
                metacluster.addEvents(clippingRightEventsToAdd)
                if discClip == True:
                    metacluster.rightClipType = 'DISC'
                else:
                    metacluster.rightClipType = 'REG'

            # If there are not discordant clippings, but metacluster has identity, return metacluster and clippings list in order to perform BLAT clippings search
            elif clippingRightEventsToAdd and not discClip and metacluster.identity:
                metaclustersWODiscClip[metacluster] = clippingRightEventsToAdd
            
            clippingLeftEventsToAdd, discClip = supportingCLIPPING(metacluster, 150, confDict, bam, normalBam, mode, 'MINUS')
            # Add clippings if there are discodant clippings.
            if discClip == True or not 'VIRUS' in confDict['targetINT2Search'] or not metacluster.identity:
                metacluster.addEvents(clippingLeftEventsToAdd)
                if discClip == True:
                    metacluster.leftClipType = 'DISC'
                else:
                    metacluster.leftClipType = 'REG'

            # If there are not discordant clippings, but metacluster has identity, return metacluster and clippings list in order to perform BLAT clippings search
            elif clippingLeftEventsToAdd and not discClip and metacluster.identity:
                if metacluster in  metaclustersWODiscClip.keys():
                    metaclustersWODiscClip[metacluster].extend(clippingLeftEventsToAdd)
                else:
                    metaclustersWODiscClip[metacluster] = clippingLeftEventsToAdd
    
    # Add clippings when there are no discordant clippings, but they have BLAT matches and clippings without blat hits but same bkp as the ones that match
    addBlatClippings(metaclustersWODiscClip, identDbPath, binId, bkpDir)

    # b. Determine bkp
    bkp(metaclusters)
    
    for metacluster in metaclusters:
        
        # c. Reconstruct bkp sequence:

        clipped_seqPlus = None
        clipped_seqMinus = None
        if metacluster.orientation == 'PLUS':
            clipped_seqPlus, clipped_seqFastaPlus = reconstructSeq(metacluster, confDict['consBkpSeq'], 'PLUS', specBkpDir)
        elif metacluster.orientation == 'MINUS':
            clipped_seqMinus, clipped_seqFastaMinus = reconstructSeq(metacluster, confDict['consBkpSeq'], 'MINUS', specBkpDir)
        elif metacluster.orientation == 'RECIPROCAL':
            clipped_seqPlus, clipped_seqFastaPlus = reconstructSeq(metacluster, confDict['consBkpSeq'], 'PLUS', specBkpDir)
            clipped_seqMinus, clipped_seqFastaMinus = reconstructSeq(metacluster, confDict['consBkpSeq'], 'MINUS', specBkpDir)
        
        # d. Determine inserted sequence bkp

        if metacluster.identity:
            if clipped_seqPlus != None:
                if not confDict['consBkpSeq'] or not clipped_seqFastaPlus:
                    # Create FASTA
                    clipped_seqPlusFastaPath = specBkpDir + '/' + binId + '_clipped_seqPlusFastaPath.fasta'
                    clipped_seqPlusFasta = formats.FASTA()
                    clipped_seqPlusFasta.seqDict['clipped_seqPlusFasta'] = clipped_seqPlus
                    clipped_seqPlusFasta.write(clipped_seqPlusFastaPath)
                elif confDict['consBkpSeq']:
                    clipped_seqPlusFastaPath = clipped_seqFastaPlus
                metacluster.intRightBkp = bkpINT(clipped_seqPlusFastaPath, identDbPath, bkpDir, metacluster.identity)
            if clipped_seqMinus != None:
                if not confDict['consBkpSeq'] or not clipped_seqFastaMinus:
                    # Create FASTA
                    clipped_seqMinusFastaPath = specBkpDir + '/' + binId + '_clipped_seqMinusFastaPath.fasta'
                    clipped_seqMinusFasta = formats.FASTA()
                    clipped_seqMinusFasta.seqDict['clipped_seqMinusFasta'] = clipped_seqMinus
                    clipped_seqMinusFasta.write(clipped_seqMinusFastaPath)
                elif confDict['consBkpSeq']:
                    clipped_seqMinusFastaPath = clipped_seqFastaMinus
                metacluster.intLeftBkp = bkpINT(clipped_seqMinusFastaPath, identDbPath, bkpDir, metacluster.identity)

    ### Do cleanup
    #TEMP
    #unix.rm([bkpDir])

def supportingCLIPPING(metacluster, buffer, confDict, bam, normalBam, mode, side):
    '''
    # Determine the area where clipping events must be searched (narrow if there are discordant clippings, wide otherwise).
    # Collect clipping events in previous area.

    Input:
        1. metacluster
        2. buffer: buffer to add in wide clippings search
        3. confDict
        4. bam
        5. normalBam
        6. mode
        7. side: orientation to look for clippings: PLUS or MINUS.
    
    Output:
        1. clippings: candidate clippings list
        2. discClip: boolean. True if there are discordant clippings, False otherwise.
    '''
    # Note: This function works but you have to allow duplicates in the clipping

    # New dictionary for performing collecting collecting clippings
    clippingConfDict = dict(confDict)
    clippingConfDict['targetEvents'] = ['CLIPPING']
    clippingConfDict['minMAPQ'] = 10
    # NOTE SR: Think and check if this is neccessary. I think it is not
    #confDict['minCLIPPINGlen'] = 2
    clippingEventsDict = {}

    ## Define region
    if side == 'PLUS':
        if metacluster.orientation != 'RECIPROCAL':
            # Determine the area where clipping events must be searched (narrow if there are discordant clippings, wide otherwise).
            binBeg, binEnd, discClip = determinePlusBkpArea(metacluster.beg, metacluster.end, metacluster.events, buffer)
        elif metacluster.orientation == 'RECIPROCAL':
            # Collect PLUS discordant events:
            reciprocalPlusEvents = [eventP for eventP in metacluster.events if eventP.orientation == 'PLUS']

            # Determine plus cluster beggining and end
            begPlus = min([eventPlus.beg for eventPlus in reciprocalPlusEvents])
            endPlus = max([eventPlus.end for eventPlus in reciprocalPlusEvents])

            # Determine the area where clipping events must be searched (narrow if there are discordatn clippings, wide otherwise).
            binBeg, binEnd, discClip = determinePlusBkpArea(begPlus, endPlus, reciprocalPlusEvents, buffer)

    elif side == 'MINUS':
        if metacluster.orientation != 'RECIPROCAL':
            # Determine the area where clipping events must be searched (narrow if there are discordatn clippings, wide otherwise).
            binBeg, binEnd, discClip = determineMinusBkpArea(metacluster.beg, metacluster.end, metacluster.events, buffer)
        elif metacluster.orientation == 'RECIPROCAL':
            # Collect MINUS discordant events
            # As some CLIPPING PLUS events could be added in previous step, now they have to be collected this way
            reciprocalMinusEvents = []
            for eventM in metacluster.events:
                if eventM.type == 'DISCORDANT':
                    if eventM.orientation == 'MINUS':
                        reciprocalMinusEvents.append(eventM)

            # Determine minus cluster beggining and end
            begMinus = min([eventMinus.beg for eventMinus in reciprocalMinusEvents])
            endPlus = max([eventMinus.end for eventMinus in reciprocalMinusEvents])

            # Determine the area where clipping events must be searched (narrow if there are discordatn clippings, wide otherwise).
            binBeg, binEnd, discClip = determineMinusBkpArea(begMinus, endPlus, reciprocalMinusEvents, buffer)

    # Collect clippings
    ref = metacluster.ref

    if mode == "SINGLE":
        clippingEventsDict = bamtools.collectSV(ref, binBeg, binEnd, bam, clippingConfDict, None, False)

    elif mode == "PAIRED":
        clippingEventsDict = bamtools.collectSV_paired(ref, binBeg, binEnd, bam, normalBam, clippingConfDict)

    # When the metacluster is RIGHT and there are some right clippings:
    if side == 'PLUS' and 'RIGHT-CLIPPING' in clippingEventsDict.keys():
        # Keep only those clipping event which have their clipping bkp is in the desired area
        clippingEventsToAdd = chooseBkpClippings(clippingEventsDict, 'RIGHT-CLIPPING', binBeg, binEnd)
        # From dictionary to list
        clippings = list(itertools.chain(*clippingEventsToAdd.values()))
        return clippings, discClip

    # When the metacluster is LEFT and there are some minus clippings:
    elif side == 'MINUS' and 'LEFT-CLIPPING' in clippingEventsDict.keys():
        # Keep only those clipping event which have their clipping bkp is in the desired area
        clippingEventsToAdd = chooseBkpClippings(clippingEventsDict, 'LEFT-CLIPPING', binBeg, binEnd)
        # From dictionary to list
        clippings = list(itertools.chain(*clippingEventsToAdd.values()))
        return clippings, discClip
    
    else:
        return [], False


def determinePlusBkpArea(beg, end, events, buffer):
    '''
    This function determines the coordinates from where clipping events should be added to a PLUS metacluster. Also, it determines if there are discordant clipping or not. 
    If there are discordant clipping, the region to look for clippings will be [the most left bkp -2 : the most left bkp +2].
    Otherwise, the the region to look for clipping events is from metacluster.beg - buffer to metacluster.end + buffer.

    Input:
        1. beg: beg position
        2. end: end position
        3. events: list of events
        4. buffer: bp of region extension when there are no discordant clippings.
    Ouput:
        1. binBeg: beg position of the region to look for clipping events.
        2. binEnd: end position of the region to look for clipping events.
        3. discClip: bool. True is there are discordatn events, False otherwise.
    '''

    plusBkp = []
    # If there are clippings in discordant events: Collect all clipping bkp genomic positions
    for discordantPlus in events:
        lastOperation, lastOperationLen = discordantPlus.cigarTuples[-1]
        # Collect all clipping bkp genomic positions
        if ((lastOperation == 4) or (lastOperation == 5)):
            plusBkp.append(discordantPlus.end)

    # The region to look for clippings will be [the most left bkp -2 : the most left bkp +2]
    if len(plusBkp) > 0:
        binBeg = min(plusBkp) - 5 if min(plusBkp) >= 5 else 0
        binEnd = max(plusBkp) + 5
        discClip = True

    
    # If there are NOT clippings in discordant events
    # The region to look for clippings will be metacluster positions + buffer
    else:
        binBeg = beg
        binEnd = end + buffer
        discClip = False

    return binBeg, binEnd, discClip

def determineMinusBkpArea(beg, end, events, buffer):
    '''
    This function determines the coordinates from where clipping events should be added to a MINUS metacluster. Also, it determines if there are discordant clipping or not. 
    If there are discordant clipping, the region to look for clippings will be [the most left bkp -2 : the most left bkp +2].
    Otherwise, the the region to look for clipping events is from metacluster.beg - buffer to metacluster.end + buffer.

    Input:
        1. beg: beg position
        2. end: end position
        3. events: list of events
        4. buffer: bp of region extension when there are no discordant clippings.

    Ouput:
        1. binBeg: beg position of the region to look for clipping events.
        2. binEnd: end position of the region to look for clipping events.
        3. discClip: bool. True is there are discordatn events, False otherwise.
    '''
    minusBkp = []
    # If there are clippings in discordant events: Collect all clipping bkp genomic positions
    for discordantMinus in events:
        firstOperation, firstOperationLen = discordantMinus.cigarTuples[0]
        # Collect all clipping bkp genomic positions
        if ((firstOperation == 4) or (firstOperationLen == 5)):
            minusBkp.append(discordantMinus.beg)

    # The region to look for clippings will be [the most left bkp -2 : the most left bkp +2]
    if len(minusBkp) > 0:
        binBeg = min(minusBkp) - 5 if min(minusBkp) >= 5 else 0
        binEnd = max(minusBkp) + 5
        discClip = True
    
    # If there are NOT clippings in discordant events
    # The region to look for clippings will be metacluster positions + buffer
    else:
        binBeg = beg - buffer if beg > buffer else 0
        # TODO SR: check as in determinePlusBkpArea
        binEnd = end
        discClip = False
    
    return binBeg, binEnd, discClip

def chooseBkpClippings(clippingEventsDict, eventType, binBeg, binEnd):
    '''
    Keep only those clipping events which clipping bkp is between input coordinates.

    Input:
        1. clippingEventsDict: clippingEventsDict[eventType] = []
        2. eventType: 'RIGHT-CLIPPING' or 'LEFT-CLIPPING'
        3. binBeg: beg position
        4. binEnd: end position
    Output:
        1. clippingEventsToAdd: same structure as input dictionary, containing only those clipping events which clipping bkp is between input coordinates.
    '''
    # eventType = 'RIGHT-CLIPPING' or 'LEFT-CLIPPING'
    clippingEventsToAdd = {}
    clippingEventsToAdd[eventType] = []
    
    ## Get clipping clusters:
    #clippingEventsDict = dict((key,value) for key, value in clippingEventsDict.items() if key == eventType)
    for clippingEvent in clippingEventsDict[eventType]:
        if eventType == 'RIGHT-CLIPPING':
            if binBeg <= clippingEvent.end <= binEnd:
                clippingEventsToAdd[eventType].append(clippingEvent)
        elif eventType == 'LEFT-CLIPPING':
            if binBeg <= clippingEvent.beg <= binEnd:
                clippingEventsToAdd[eventType].append(clippingEvent)

    return clippingEventsToAdd

def addBlatClippings(metaclustersWODiscClip, db, binId, outDir):
    '''
    Perform BLAT search: add clipping events with BLAT hits and clipping events sharing their position, even they have no matches.
    Input:
        1. metaclustersWODiscClip: dictionary -> metaclustersWODiscClip[metacluster] = [clippings]. key -> metacluster object; value -> list of candidate clipping events 
        2. db: Path to fasta file to be used as blat reference.
        3. binId
        4. outDir
    Output:
        Doesn't return anything.
        Add clippings to metacluster if they fulfill the above conditions.
    '''
    # Write fasta with events collected in a wide region
    clippingEventsToAdd = list(itertools.chain(*metaclustersWODiscClip.values()))
    clippingsFasta = writeClippingsFasta(clippingEventsToAdd, binId, outDir)
    # Align clipped sequences with BLAT against db
    blatArgs = {}
    blatArgs['tileSize'] = 7
    outName = binId + '_clippingsBlat'
    pslPath = alignment.alignment_blat(clippingsFasta, db, blatArgs, outName, outDir)
    # Make dictionary from blat results.
    pslDict = formats.pslQueryRefDict(pslPath)
    matchClippings = collectMatchClippings(metaclustersWODiscClip, pslDict)
    if matchClippings != {}:
        # Collect those clippings that have their bkp in same position as ones in BLAT.
        clippings2Add = collectClipBkpMatch(matchClippings, clippingEventsToAdd)
        # Add all clippings to metacluster
        if clippings2Add:
            for metacluster, clippings in clippings2Add.items():
                metacluster.addEvents(clippings)
                if metacluster.orientation == 'PLUS':
                    metacluster.rightClipType = 'BLAT'
                elif metacluster.orientation == 'MINUS':
                    metacluster.leftClipType = 'BLAT'

def writeClippingsFasta(clippings, binId, outDir):
    '''
    Perform a BLAT search with clipping events, select bkp area of those that match and pick all clipping event which bkp is in this area.

    Input:
        1. clippingEventsDict: clippingEventsDict[eventType] = []
        2. eventType: 'RIGHT-CLIPPING' or 'LEFT-CLIPPING'
        3. identity: metacluster identity
        4. ID: metacluster ID
        5. db: reference DB
        6. outDir: output directory
    '''

    ## 1. Generate fasta containing soft clipped sequences
    clippedFasta = events.collect_clipped_seqs(clippings)

    ## 2. Write clipped sequences into fasta
    filePath = outDir + '/' + binId + '_clippings.fasta'
    clippedFasta.write(filePath, 'append', False)

    return filePath


def collectMatchClippings(metaclustersWODiscClip, pslDict):
    '''
    Collect those clippings with hits in BLAT search whose hit match with metacluster identities.

    Input:
        1. metaclustersWODiscClip: dictionary -> metaclustersWODiscClip[metacluster] = [clippings]. key -> metacluster object; value -> list of candidate clipping events
        2. pslDict: dictionary -> pslDict[qName] = tName 
    
    Output:
        This function changes clippingEvent.blatIdentity attribute: True if event has blat hits that match with metacluster identity, False otherwise.
        1. matchClippings. Dictionary with same structure as input one, but containing as values a list with only those clippings that have blat match
    '''
    matchClippings = {}
    for metaclusterWODiscClip, clippings in metaclustersWODiscClip.items():
        for clip in clippings:
            if clip.readName in pslDict.keys():
                # If metacluster.identity is a string
                # Check if metacluster.identity is in blat hits list. Identity has partial name (i.e Hepadnaviridae) wheter pslDict complete names (i.e. Hepadnaviridae|KR811803.1).
                if type(metaclusterWODiscClip.identity) is str:
                    if any(metaclusterWODiscClip.identity in iden for iden in pslDict[clip.readName]):
                        if metaclusterWODiscClip in matchClippings.keys():
                            matchClippings[metaclusterWODiscClip].append(clip)
                            clip.blatIdentity = True
                        else:
                            matchClippings[metaclusterWODiscClip] = []
                            matchClippings[metaclusterWODiscClip].append(clip)
                            clip.blatIdentity = True
                # If metacluster.identity is a list
                # Check if any element of metacluster.identity is in blat hits list. Identity has partial names (i.e Hepadnaviridae) wheter pslDict complete names (i.e. Hepadnaviridae|KR811803.1).
                elif type(metaclusterWODiscClip.identity) is list:
                    if [i for e in metaclusterWODiscClip.identity for i in pslDict[clip.readName] if e in i]:
                        if metaclusterWODiscClip in matchClippings.keys():
                            matchClippings[metaclusterWODiscClip].append(clip)
                            clip.blatIdentity = True
                        else:
                            matchClippings[metaclusterWODiscClip] = []
                            matchClippings[metaclusterWODiscClip].append(clip)
                            clip.blatIdentity = True                     

    return matchClippings
                
def collectClipBkpMatch(matchClippings, clippingEventsToAdd):
    '''
    From matchClippings that contains clipping events with blat hits and clippingEventsToAdd list containing all candidate clipping events (including previous ones)
    Make a new dictionary containing all clipping from matchClippings + clippings from clippingEventsToAdd with same bkp as those in matchClippings

    Input:
        1. matchClippings: Dictionary with same structure as input one, but containing as values a list with only those clippings that have blat match
        2. clippingEventsToAdd: clipping events list
    
    Output:
        1. clippings2Add: matchClippings + clipping match bkp events.
    '''
    
    # Collect bkp position of BLAT hits clippping events
    clippings2Add = {}
    for metacluster, clippings in matchClippings.items():
        bkp = []
        for matchClipping in clippings:
            bkp.append(matchClipping.readBkp)

        # Make coordinates of the region where bkp should be
        binBeg = min(bkp) - 5 if min(bkp) >= 5 else 0
        binEnd = max(bkp) + 5

        # Collect clipping events whose bkp is in desired region
        for clipping in clippingEventsToAdd:
            if clipping.readBkp >= binBeg and clipping.readBkp <= binEnd and clipping in clippings:
                if metacluster in clippings2Add.keys():
                    clippings2Add[metacluster].append(clipping)
                    
                else:
                    clippings2Add[metacluster] = []
                    clippings2Add[metacluster].append(clipping)
    
    return clippings2Add


def reconstructSeq(metacluster, consSeq, orientation, outDir):
    '''
    Reconstruct bkp sequence.

    Input:
        1. metacluster
        2. consSeq: boolean. If True make consensus sequence. Otherwise, choose a representative alignment as bkp sequence reconstruction.
        3. orientation: 'PLUS' or 'MINUS'
        4. outDir
    
    Output:
        Fill following metacluster attributes: metacluster.rightSeq and metacluster.leftSeq
        1. clipped_seq: Path to consensus fasta file if it was made. Otherwise, clipped part of representative sequence.
        2. consFastaBool: boolean. True is consensus fasta file was made, False otherwise.
    '''
    # Collect discordant clipping events in dictionary, having clipping length as value
    discClip = {}
    clippingsDisc = {}
    clipped_seqFasta = None
    for event in metacluster.events:
        if event.type == 'DISCORDANT':
            if orientation == 'PLUS':
                lastOperation, lastOperationLen = event.cigarTuples[-1]
                # Collect all clipping bkp genomic positions
                if ((lastOperation == 4) or (lastOperation == 5)):
                    discClip[event.readName] = lastOperationLen
            elif orientation == 'MINUS':
                firstOperation, firstOperationLen = event.cigarTuples[0]
                # Collect all clipping bkp genomic positions
                if ((firstOperation == 4) or (firstOperationLen == 5)):
                    discClip[event.readName] = firstOperationLen
    # If there are discordant clipping events      
    if discClip:
        # Collect clipping objects that are the same alignment as discordant object:
        clippingsDisc = {}
        for eventC in metacluster.events:
            if eventC.type == 'CLIPPING' and eventC.readName in discClip.keys():
                clippingsDisc[eventC] = discClip[eventC.readName]

    
    if clippingsDisc: # This can be empty if clipping of discordant is too small.
        if consSeq and len(clippingsDisc) > 1 and len(clippingsDisc) < 1000: # Make consensus sequence with discordant clipping events. len(clippingsDisc) > 1 in order to avoid long lasting muscle runs
            clipped_seq, clipped_seqFasta = conSeq(metacluster, clippingsDisc, orientation, outDir)
        else: # Make representative sequence with discordant clipping events
            clipped_seq = repreSeq(metacluster, orientation, clippingsDisc)
        return clipped_seq, clipped_seqFasta


    # If there are no discordant clipping events  
    else:
        # Collect those clippings with blat hits in dictionary, having clipping length as value
        clippingsBlat = {}
        for clippingB in metacluster.events:
            if clippingB.type == 'CLIPPING':
                if (orientation == 'PLUS' and clippingB.clippedSide == 'right') or (orientation == 'MINUS' and clippingB.clippedSide == 'left'):
                    if clippingB.blatIdentity == True:
                        clippingsBlat[clippingB] = clippingB.cigarTuples[-1][1]

        # If there are clippings with BLAT hits
        if clippingsBlat:
            if consSeq and len(clippingsBlat) > 1 and len(clippingsBlat) < 1000: # Make consensus sequence with BLAT clipping events. len(clippingsDisc) > 1 in order to avoid long lasting muscle runs
                clipped_seq, clipped_seqFasta = conSeq(metacluster, clippingsBlat, orientation, outDir)
            else: # Make representative sequence with BLAT clipping events
                clipped_seq = repreSeq(metacluster, orientation, clippingsBlat)
            return clipped_seq, clipped_seqFasta

        # If there are no clippings with BLAT hits
        else:
            # Collect metacluster clippings
            clippings = {}
            for clipping in metacluster.events:
                if clipping.type == 'CLIPPING':
                    if (orientation == 'PLUS' and clipping.clippedSide == 'right') or (orientation == 'MINUS' and clipping.clippedSide == 'left'):
                        clippings[clipping] = clipping.cigarTuples[-1][1]
            # If there are clippings
            if clippings:
                # NOTE SR: this sequence will be less relayable
                if consSeq and len(clippings) > 1 and len(clippings) < 1000: # Make consensus sequence with clipping events. len(clippingsDisc) > 1 in order to avoid long lasting muscle runs
                    clipped_seq, clipped_seqFasta = conSeq(metacluster, clippings, orientation, outDir)
                else: # Make representative sequence with clipping events
                    clipped_seq = repreSeq(metacluster, orientation, clippings)
                return clipped_seq, clipped_seqFasta

            # If there are no clippings, there are not representative sequence.
            else:
                return None, None

def repreSeq(metacluster, orientation, clippings):
    '''
    Choose an alignment as representative bkp reconstruction.
    Criteria: Choose the one with longest clipping sequence.

    Input:
        1. metacluster:
        2. orientation: bkp orientation -> 'PLUS' or 'MINUS'
        3. clippings: dictionary -> key: clipping events that are candidates of be the representative alignment
    
    Output:
        This function fill metacluster.rightSeq when orientation == 'PLUS' and metacluster.leftSeq when orientation == 'MINUS'
        1. largestClipping.clipped_seq(): clipping part of representative sequence.
    '''
    
    # Choose the one with maximum clipping lenght
    largestClipping = max(clippings.items(), key=operator.itemgetter(1))[0]
    # Make the representative sequence
    if orientation == 'PLUS':
        metacluster.repreRightSeq = largestClipping.ref_seq() + '[INT]>' + largestClipping.clipped_seq()
    elif orientation == 'MINUS':
        metacluster.repreLeftSeq = largestClipping.clipped_seq() + '<[INT]' + largestClipping.ref_seq()
    return largestClipping.clipped_seq()

def conSeq(metacluster, clippings, orientation, outDir):
    '''
    Perform consensus sequence from clipping events

    Input:
        1. metacluster
        2. clippings: dictionary -> key: clipping events that are candidates of be the representative alignment
        3. orientation: bkp orientation -> 'PLUS' or 'MINUS'
        4. outDir
    
    Output:
        Fill metacluster attributes: metacluster.rightSeq and metacluster.leftSeq
        1. intConsensusSeq: Insertion part of consensus sequence
        2. intConsensusPath: Path to consensus fasta file
    '''
    refConsensusSeq = makeConsSeqs([*clippings], 'REF', outDir)[1]
    intConsensusPath, intConsensusSeq = makeConsSeqs([*clippings], 'INT', outDir)
    if intConsensusSeq != None:
        if orientation == 'PLUS':
            metacluster.consRightSeq = refConsensusSeq + '[INT]>' + intConsensusSeq
        elif orientation == 'MINUS':
            metacluster.consLeftSeq = intConsensusSeq + '<[INT]' + refConsensusSeq
    return intConsensusSeq, intConsensusPath

def makeConsSeqs(clippingEvents, seqSide, outDir):
    '''
    Make consesus sequence of one of the sides of the clipping cluster.

    Input:
        1. CLIPPING_cluster
        2. clippedSide: 'left' or 'right'
        3. seqSide: 'INT' if the consensus of the integrated sequence is wanted or 'REF' if the consensus of the reference is wanted
    
    Output:
        1. consensusPath: consensus file
        2. consensusSeq: consensus sequence
    '''
    consensusPath = None
    consensusSeq = None

    #clippingEvents = [event for event in CLIPPING_cluster.events if event.clippedSide == clippedSide]

    if len (clippingEvents) > 0:

        consensusPath, consensusSeq = clippingConsensusSeq(clippingEvents, clippingEvents[0].id, seqSide, outDir)
    
    return consensusPath, consensusSeq


def clippingConsensusSeq(clippingEvents, CLIPPING_clusterID, seqSide, outDir):

    # Retrieve fasta file with sequence from match or clipped side of clipping reads
    supportingReadsFasta = clippingSeq(clippingEvents, CLIPPING_clusterID, seqSide, outDir)

    # Consensus from the previous fasta
    consensusPath, consensusSeq = assembly.getConsensusSeq(supportingReadsFasta, outDir)

    return consensusPath, consensusSeq


def clippingSeq(clippingEvents, CLIPPING_clusterID, seqSide, outDir):
    '''
    Retrieve fasta file with sequence from match or clipped side of clipping reads
    '''

    fastaObj = formats.FASTA()
    fastaDict = {}
                
    # Determine bkp
    for event in clippingEvents:
        
        # Si queremos sacar la secuencia del lado de la integracion:
        if seqSide == 'INT':
            fastaDict[event.readName] = event.clipped_seq()

        elif seqSide == 'REF':
            fastaDict[event.readName] = event.ref_seq()
        
        fastaObj.seqDict = fastaDict

    fastaPath = outDir + '/' + str(CLIPPING_clusterID) +'_'+ str(seqSide) +'_supportingReads.fa'
    fastaObj.write(fastaPath)

    return fastaPath


def bkpINT(fastaPath, db, outDir, identity):
    '''
    Align sequence with minimap2 and get bkp pos in reference. Keep several hits, but only those ones that match with metacluster identity.

    Input:
        1. fastaPath: Path to FASTA file
        2. db: Path to reference 
        3. outDir
        4. identity: metacluster indentity
    
    Output:
        1. intBkp -> dictionary: intBkp[referenceName] = referencePosition
    '''

    #indexDbSpecificIdentity = databases.buildIdentityDb(metacluster, db, outDir)   

    PAF_file = sequences.getPAFAlign(fastaPath, db, outDir)
    PAFObj = formats.PAF()
    PAFObj.read(PAF_file)

    intBkp = {}
    if not os.stat(PAF_file).st_size == 0:
        for alig in PAFObj.alignments:
            # Check identity
            if alig.tName.split('|')[0] in identity:
                intBkp[alig.tName.split('|')[1]] = alig.tBeg
        #intBkp = [line.tBeg for line in PAFObj.alignments][0]
    else:
        intBkp = {}

    return intBkp
