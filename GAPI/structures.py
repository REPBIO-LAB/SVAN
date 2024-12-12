'''
Module 'structures' - Contains functions and classes to organize data into more complex data structures
'''

## DEPENDENCIES ##
# External
import multiprocessing as mp

# Internal
from GAPI import log
from GAPI import gRanges

## FUNCTIONS ##

def create_bin_database(refLengths, eventsDict):
    '''
    Organize genome wide events into a set of bin databases, one per reference

    Input:
        1. refLengths: dictionary containing references as keys and their lengths as values
        2. eventsDict: nested dictionary containing:
            * FIRST LEVEL KEYS:
                - REF_1
                - ...

                * SECOND LEVEL KEYS:
                    - EVENT_TYPE_1 -> list of objects
                    - ...        
    Output:
        1. wgBinDb: dictionary containing references as keys and the corresponding 'bin_database' as value
    '''    
    # Initialize dict
    wgBinDb = {}

    # For each reference
    for ref, refLen in refLengths.items():
        
        # Skip if no events in that particular ref
        if ref not in eventsDict:
            continue

        # Define bin sizes
        #binSizes = [100, 1000, 10000, 100000, int(refLen)]
        #binSizes = [1000000, 10000000, int(refLen)]
        binSizes = [10000, 100000, 1000000, int(refLen)]
        
        # Create bin database for reference
        binDb = create_bin_database_interval(ref, 0, refLen, eventsDict[ref], binSizes)      

        # Add reference bin database to the dictionary
        wgBinDb[ref] = binDb
    
    return wgBinDb


def create_bin_database_parallel(refLengths, eventsDict, threads):
    '''
    Organize genome wide events into a set of bin databases, one per reference

    Input:
        1. refLengths: dictionary containing references as keys and their lengths as values
        2. eventsDict: nested dictionary containing:
            * FIRST LEVEL KEYS:
                - REF_1
                - ...

                * SECOND LEVEL KEYS:
                    - EVENT_TYPE_1 -> list of objects
                    - ...
        3. threads: number of threads used to parallelize the bin database creation
        
    Output:
        1. wgBinDb: dictionary containing references as keys and the corresponding 'bin_database' as value
    '''    
    ## 1. Create list of tuples for multiprocessing
    # Each tuple will contain all the variables needed for loading the events 
    # into a the whole genome bin database
    tupleList = []

    # For each reference
    for ref, refLen in refLengths.items():
        
        # Skip if no events in that particular ref
        if ref not in eventsDict:
            continue

        # Define bin sizes
        binSizes = [10000, 100000, 1000000, int(refLen)]

        # Add to the list of tuples
        fields = (ref, 0, refLen, eventsDict[ref], binSizes)
        tupleList.append(fields)
    
    ## 2. Create bin database per chromosome
    pool = mp.Pool(processes=threads)
    databases = pool.starmap(create_bin_database_interval, tupleList)
    pool.close()
    pool.join()

    ## 3. Organize bin databases into a dictionary
    wgBinDb = {}

    for binDb in databases:
        wgBinDb[binDb.ref] = binDb
    
    return wgBinDb

def create_bin_database_interval(ref, beg, end, eventsDict, binSizes):
    '''
    Organize events into a bin database. Events are removed from the original dict once they are incorporated into the bin database

    Input:
        1. ref: reference/chromosome
        2. beg: bin begin coordinate
        3. end: bin end coordinate
        4. eventsDict:  Dictionary containing:

            - EVENT_TYPE_1 -> list of objects
            - ...
    
        5. binSizes: list of bin sizes that will be used to create a bin database 

    Output:
        1. binDb: 'bin_database' instance containing all the input events organized in genomic bins
    '''            
    ## Initiate bin database
    binDb = bin_database(ref, beg, end, binSizes)

    ## For each type of input event
    for eventType, events in eventsDict.items():

        # Add all the events from the given event type to the bin database
        binDb.add(events, eventType)
    
    return binDb

## CLASSES ##
class bin_database():
    '''
    Database to organize a set of events into a hierarchy of genomic bins
    '''
    def __init__(self, ref, beg, end, binSizes):
        '''
        Organize events into a hierarchy of genomic bins

        Input:
            1. binSizes: list of bin sizes 
        '''
        self.ref = ref
        self.beg = int(beg)  
        self.end = int(end)
        self.binSizes = sorted(binSizes)
        self.eventTypes = []
        self.data = {}
        
        ## Initialize one bin dictionary per size
        for binSize in self.binSizes:
            self.data[binSize] = {}

    def add(self, events, eventType):
        '''
        Add events into a hierarchy of genomic bins

        Input:
            1. events: list of objects. Every object must have 'beg' and 'end' attributes
            2. eventType: type of events (DEL, INS, CLIPPING, ...)
        '''        
        ## 1. Add event type to event type list
        if eventType not in self.eventTypes:
            self.eventTypes.append(eventType)

        ## 2. Allocate each event into a genomic bin 
        # For each event
        for event in events:
        
            # For each bin size (from smaller to bigger bin sizes)
            for binSize in self.binSizes:

                # Determine to what bin index the event belongs
                binIndexBeg = int(event.beg / binSize)
                binIndexEnd = int(event.end / binSize)

                # A) Event fits in one bin
                if (binIndexBeg == binIndexEnd):

                    # a) First event into that bin -> Initialize dict and bin  
                    if binIndexBeg not in self.data[binSize]:
                        self.data[binSize][binIndexBeg] = {}
                        self.data[binSize][binIndexBeg][eventType] = events_bin([event])                

                    # b) First event of that type in this bin -> Initialize bin
                    elif eventType not in self.data[binSize][binIndexBeg]:
                        self.data[binSize][binIndexBeg][eventType] = events_bin([event]) 

                    # c) There are already events of this type in this bin -> Add event to the bin
                    else:
                        self.data[binSize][binIndexBeg][eventType].add([event])
       
                    # Do not check other bin sizes once event allocated in a bin
                    break

                ## B) Event spans several bins. Try with the next bin size 
                else:
                    continue

        ## 3. For each bin sort the events in increasing coordinates order
        for binSize in self.data.keys():
            for binIndex in self.data[binSize].keys():
                if eventType in self.data[binSize][binIndex]:
                    self.data[binSize][binIndex][eventType].sort()
    
    def remove(self, events, eventType):
        '''
        Remove events from a hierarchy of genomic bins

        Input:
            1. events: list of objects
            2. eventType: type of events (DEL, INS, CLIPPING, ...)
        '''        

        # For each event
        for event in events:
        
            # For each bin size (from smaller to bigger bin sizes)
            for binSize in self.binSizes:
    
                # Determine to what bin index the event belongs
                binIndexBeg = int(event.beg / binSize)
                binIndexEnd = int(event.end / binSize)

                ## A) Target bin found if:
                # - Event fits in a single bin using current bin size
                # - Bin index available AND
                # - There are events of the target event type in the target bin 
                if (binIndexBeg == binIndexEnd) and (binIndexBeg in self.data[binSize]) and (eventType in self.data[binSize][binIndexBeg]):

                    binObj = self.data[binSize][binIndexBeg][eventType]
                    
                    # Remove event from bin
                    binObj.remove([event])

                    # Do not check other bin sizes once event has been removed
                    break

                ## B) Event does not fit. Try with the next bin size 
                else:
                    continue

    def collect(self, eventTypes):
        '''
        Collect all the events of target event types that are stored 
        in the bin database structure
        
         Input:
            1. eventTypes: list containing target event types

         Output:
            2. events. List of events
        '''  
        events = []

        # For each bin size
        for binSize in self.binSizes:

            # For each bin
            for binIndex in self.data[binSize].keys():
                
                # For each target event type
                for eventType in eventTypes:

                    # There are events of the target event type in the bin 
                    if eventType in self.data[binSize][binIndex]:

                        # Add events to the list
                        events = events + self.data[binSize][binIndex][eventType].events

        ## Sort events by begin coordinate
        events.sort(key=lambda event: event.beg)

        return events

    def collect_bin(self, binSize, binIndex, eventTypes):
        '''
        Collect all the events of target event types that are stored 
        in a particular bin 

         Input:
            1. binSize: bin size corresponding to the target bin index
            2. binIndex: target bin index
            3. eventTypes: list containing target event types

         Output:
            2. events. List of events
        '''  
        events = []

        ## Check if bin database contains target bin
        if (binSize in self.data) and (binIndex in self.data[binSize]):
    
            # For each target event type
            for eventType in eventTypes:

                # There are events of the target event type in the bin 
                if eventType in self.data[binSize][binIndex]:

                    # Add events to the list
                    events = events + self.data[binSize][binIndex][eventType].events

        ## Sort events by begin coordinate
        events.sort(key=lambda event: event.beg)

        return events

    def collect_interval(self, beg, end, eventTypes):
        '''
        Collect all the events of target event types that are located with the input interval

         Input:
            1. beg: interval begin
            2. end: interval end
            3. eventTypes: list containing target event types. If 'ALL' all event types stored in the target bins will be retrieved

         Output:
            1. events. List of lists. Each list corresponds to one overlapping event and is composed by 4 elements: 
                1. Overlapping event
                2. Number of overlapping base pairs
                3. Percentage of base pairs of the input interval that are overlapping  
                4. Tuple with input interval coordinates overlapping with the event
        '''  
        ## 1. Determine genomic bins spanned by the input interval
        minBinSize = self.binSizes[0]
        indexStart = int(beg / minBinSize)
        indexEnd = int(end / minBinSize)
            
        ## 2. Collect all the events potentially within input interval 
        allEvents = []

        # For each bin spanning the interval
        for index in range(indexStart, indexEnd + 1, 1):

            # Collect events at target bin and overlapping bins at higher levels of the database
            events = self.traverse(index, minBinSize, eventTypes)
            allEvents = allEvents + events
            
        # Remove redundant events
        allEvents = list(set(allEvents))

        ## 3. Make list of events within input interval
        events = []

        # For each event
        for event in allEvents:

            # Assess if located within the interval
            overlap, overlapLen, coord = gRanges.overlap_extended(beg, end, event.beg, event.end)

            # Compute percentage of overlap
            intervalLen = end - beg + 1
            overlapPerc = float(overlapLen) / intervalLen * 100
            
            # Overlap found -> Add to the list the tuple
            if overlap:
                events.append([event, overlapLen, overlapPerc, coord])
        
        return events
    
    def traverse(self, rootIndex, rootSize, eventTypes):
        '''
        Traverse bin structure starting in a root bin and going through all the bins located at upper levels
        in the hierarchy. Collect events from all the visited bins. E.g:

        <-----------------1---------------->
        <-------2--------><-------3-------->
        <---4---><---5---><---6---><---7---> 
        # BinIndex=4; Output: events in bins (4, 2, 1)
        # BinIndex=7; Output: events in bins (7, 3, 1)
        # BinIndex=2; Output: events in bins (2, 1)

        Input:
            1. rootIndex: root bin index
            2. rootSize: window size/level where the root index is located 
            3. eventTypes: list containing target event types. If 'ALL' all event types stored in the target bins will be retrieved

        Output:
            1. events: list of events 
        '''    
        # 'eventTypes' set as 'ALL' -> All distinct event types stored in the bin database will be retrieved
        if eventTypes == 'ALL':
            eventTypes = self.eventTypes

        ### Initialize events list adding the events from the root bin
        try:
            events = self.collect_bin(rootSize, rootIndex, eventTypes)
            
            ### Select upper windows sizes/levels
            upperSizes = [ binSize for binSize in self.binSizes if binSize > rootSize]
            
            # For each upper window size
            for upperSize in upperSizes:
                    
                # Compute corresponding upper bin index
                upperIndex = int(rootIndex * rootSize / upperSize)

                # Collect events in the upper bin
                upperBinEvents = self.collect_bin(upperSize, upperIndex, eventTypes)

                ## Add upper bin´s events to the list
                events = events + upperBinEvents
            
        except KeyError:
            events =  []

        return events

    def nbEvents(self):
        '''
        Compute the number of events composing the bin dictionary structure
        '''
        totalNbEvents = 0
        nbEventsBinSizes = {}

        for binSize in self.binSizes:
            nbEventsBinSizes[binSize] = 0

            for binIndex in self.data[binSize].keys():
                for eventType in self.data[binSize][binIndex].keys():
                    totalNbEvents += self.data[binSize][binIndex][eventType].nbEvents()
                    nbEventsBinSizes[binSize] += self.data[binSize][binIndex][eventType].nbEvents()

        return totalNbEvents, nbEventsBinSizes


class events_bin():
    '''
    Contain a set of events of the same type in a genomic bin
    '''
    def __init__(self, events):
        self.events = events

    def add(self, events):
        '''
        Add input events to the bin´s list of events

        Input:
            1. events: List of events to add to the bin 
        '''
        self.events = self.events + events

    def remove(self, events):
        '''
        Remove input events from the bin´s list of events

        Input:
            1. events: List of events to remove from the bin 
        '''

        ## Make list with event id´s to remove
        ids = []

        for event in events:

            ids.append(event.id)

        ## Remove events with matching ids from the bin
        self.events = [ event for event in self.events if event.id not in ids ]
            
    def nbEvents(self):
        '''
        Compute the number of events composing the bin
        '''
        return len(self.events)

    def sort(self):
        '''
        Sort events in increasing coordinates order
        '''
        self.events.sort(key=lambda event: event.beg)

def merge_dictionaries(dictionaries):
    '''
    Takes as input a list of dictionaries and merge them by their keys into a single dictionary. 
    Dictionary values must be lists
    '''
    ## 1. Generate set containing all possible dictionary keys
    keys = set().union(*(dictionary.keys() for dictionary in dictionaries))

    ## 2. Initialize empty dictionary
    outDict = {} 

    for key in keys:
        outDict[key] = []   

    ## 3. Fill dictionary 
    # For each input dict
    for dictionary in dictionaries:

        # For each key and list
        for key, elements in dictionary.items():
            outDict[key] = outDict[key] + elements

    return outDict

def dict2list(dictionary):
    '''
    Takes as input a dictionary containing lists as values and return a single list containing all 
    the elements
    '''
    allItems = []

    for items in dictionary.values():
        allItems = allItems + items 

    return allItems