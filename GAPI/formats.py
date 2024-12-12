'''
Module 'formats' - Contains classes for dealing with file formats such as fasta, bed, vcf, etc...
'''

## DEPENDENCIES ##
# External
import itertools
import sys
import time
import re
import multiprocessing as mp

# Internal
from GAPI import log
from GAPI import gRanges
from GAPI import structures

## FUNCTIONS ##
def chrom_lengths_index(index):
    '''
    Read FASTA index and build a dictionary containing chromosome lengths

    Input:
        1. index: FASTA file index generated with samtools faidx
    Output:
        1. chromLen: Dictionary containing the length for each chromosome 
    '''

    chromLengths = {}

    with open(index) as indexFile:
        for line in indexFile:
            line = line.rstrip().split("\t")
            ref = line[0]
            length = line[1]
            chromLengths[ref] = int(length)

    return chromLengths

def merge_FASTA(FASTA_list):
    '''
    Merge a list of FASTA objects into a single one

    Input:
        1. FASTA_list: List of FASTA objects to be merged
    Output:
        1. FASTA_merged: FASTA object resulting from merging
    '''
    ## Initialize output FASTA
    FASTA_merged = FASTA()

    ## For each FASTA object
    for FASTA_obj in FASTA_list:

        # For each sequence in the FASTA
        for seqId, seq in FASTA_obj.seqDict.items():
                    
            # Add sequence to the merged FASTA
            FASTA_merged.seqDict[seqId] = seq 

    return FASTA_merged

def bed2binDb(bedPath, refLengths, threads):
    '''
    Organize features in a bed file into a whole genome bin database. 
    
    Bed file must contain at least 4 fields: ref, beg, end and name (feature name). Extra fields
    will not be considered

    Input:
        1. bedPath: path to bed file
        2. refLengths: Dictionary containing reference ids as keys and as values the length for each reference
        3. threads: number of threads used to parallelize the bin database creation

    Output:
        1. wgBinDb: dictionary containing references as keys and the corresponding 'bin_database' as value
    '''    
    ## Read bed
    bed = BED()
    targetRefs = list(refLengths.keys())
    bed.read(bedPath, 'nestedDict', targetRefs)
    
    ## Create bin database
    wgBinDb = structures.create_bin_database(refLengths, bed.lines)

    return wgBinDb

def INS2binDb(VCFs, refLengths, threads):
    '''
    Organize INS events from a set of input VCFs into a bin database
    
    Input:
        1. VCFs: list of VCF of objects containing INS events
        2. refLengths: Dictionary containing reference ids as keys and as values the length for each reference
        3. threads: number of threads used to parallelize the bin database creation

    Output:
        1. wgBinDb: dictionary containing references as keys and the corresponding 'bin_database' as value
    '''
    ## For each insertion event compute beg and end coordinates:
    for VCF in VCFs:

        for variant in VCF.variants:

            variant.beg, variant.end = variant.pos_interval()

    ## Organize INS events into a dict
    insDict = INS2Dict(VCFs)

    ## Create bin database
    wgBinDb = structures.create_bin_database_parallel(refLengths, insDict, threads)

    return wgBinDb

def INS2Dict(VCFs):
    '''
    Organize INS events from a set of input VCFs into a dictionary 

    Input:
        1. VCFs: List of VCF objects

    Output:
        1. VCFs_dict: Nested dictionary with first level keys as chromosomes, second level as insertion type and list of INS variants as values 
    '''
    outDict = {}

    for VCF in VCFs:

        for variant in VCF.variants:
            ref = variant.chrom
            iType = variant.info['ITYPE'] + '-' + variant.info['FAM'] if 'FAM' in variant.info else variant.info['ITYPE'] 

            if ref not in outDict:
                outDict[ref] = {}
                outDict[ref][iType] = [variant]
                continue
            
            if iType not in outDict[ref]:
                outDict[ref][iType] = [variant]
                continue

            outDict[ref][iType].append(variant)

    return outDict

## CLASSES ##
class FASTA():
    '''
    Class for dealing with files in FASTA format
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.seqDict = {}

    def read(self, filePath):
        '''
        FASTA file reader. Read and add data to a dictionary with the following format:
            - Key. Sequence identifier
            - Value. Sequence 
        '''
        fastaFile = open(filePath)

        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in itertools.groupby(fastaFile, lambda line: line[0] == '>'))

        for header in faiter:
            # drop the >
            header = header.__next__()[1:].strip()
            
            # drop the info
            header = header.split(' ')[0]

            # join all sequence lines to one.
            seq = ''.join(s.strip() for s in next(faiter))
            self.seqDict[header] = seq

    def write(self, filePath):
        '''
        FASTA file writer. Write data stored in the dictionary into a FASTA file
        Mode: write -> write new file. append -> append to existing file or create if tit doesnt exist.
        '''

        fastaFile = open(filePath, 'w')

        for header, seq in self.seqDict.items():
            header = '>' + header

            fastaFile.write("%s\n" % header)
            fastaFile.write("%s\n" % seq)

        # Close output fasta file
        fastaFile.close()

    def retrieve_seqs(self, targetNames):
        '''
        Retrieve set of sequences from fasta file

        Input:
            1. targetNames: list of read ids to be retrieved
        
        Output:
            2. outDict: dictionary containing sequences
        '''
        outDict = {readName: self.seqDict[readName] for readName in targetNames if readName in self.seqDict}
        
        if safetyLock:
            callers.lock.release()

    def retrieve_seqs(self, targetNames):
        '''
        Retrieve set of sequences from fasta file

        Input:
            1. targetNames: list of read ids to be retrieved
        
        Output:
            2. outDict: dictionary containing sequences
        '''
        outDict = {readName: self.seqDict[readName] for readName in targetNames if readName in self.seqDict}
        
        return outDict

class FASTQ():
    '''
    Class for dealing with files in FASTQ format
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.seqDict = {}

    def write(self, filePath):
        '''
        FASTQ file writer. Write data stored in the dictionary into a FASTA file
        '''
        # Open fastq file
        fastqFile = open(filePath, 'w')

        # Write FASTQ lines
        for FASTQ_entry in self.seqDict.values(): 

            seqId = '@' + FASTQ_entry.seqId
            description = '+' + FASTQ_entry.description

            fastqFile.write("%s\n" % seqId)
            fastqFile.write("%s\n" % FASTQ_entry.seq)
            fastqFile.write("%s\n" % description)
            fastqFile.write("%s\n" % FASTQ_entry.qual)

        # Close output fasta file
        fastqFile.close()

    def add(self, fastqLine):
        '''
        Add sequence to the fastq object
        '''
        self.seqDict[fastqLine.seqId] = fastqLine
        

class FASTQ_entry():
    '''
    FASTQ entry class 
    '''

    def __init__(self, seqId, seq, description, qual):
        '''
        Initialize fastq line
        '''
        self.seqId = seqId
        self.seq = seq
        self.description = description
        self.qual = qual

class BED():
    '''
    Class for dealing with files in BED format. 
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.lines = []
        self.structure =  None

    def read(self, filePath, structure, targetRefs):
        '''
        BED file reader. Read and store bed lines into a data structure

        Input:
            1) filePath: path to bed file
            2) structure: data structure where BED lines will be stored. 3 Possible structures:
                - 'List': lines saved in a list
                - 'Dict': lines saved in a dictionary where each key will correspond to a reference and the corresponding value will be the list of lines in that reference
                - 'nestedDict': lines saved in a nested dictionary where first level keys correspond to the references and second level keys to bed entry names

            3) targetRefs: list containing target references. Entries with refs not included in the list will not be loaded. 
                            If 'None' all the entries will be loaded
        Initialize lines attribute as output
        '''
        self.structure = structure

        # a) Organize BED entries into a list
        if (self.structure == 'List'):
            self.lines = self.organize_list(filePath, targetRefs)

        # b) Organize BED entries into a dict
        elif (self.structure == 'Dict'):
            self.lines = self.organize_dict(filePath, targetRefs)

        # c) Organize BED entries into a nested dict
        elif (self.structure == 'nestedDict'):
            self.lines = self.organize_nestedDict(filePath, targetRefs)

        # d) Unkown data type structure provided
        else:
            print('[ERROR] Bed file reader. Unknown structure provided: ', structure)
            sys.exit(1)

    def write(self, outPath, optional):
        '''
        BED file writer. Write bed entries into bed file

        Input:
            1) outPath: path to output bed file
            2) optional: list with optional fields to include into the output file. Empty if no optional fields

        '''
        ## Collect all the entries into a list
        # a) Entries organized into a list
        if (self.structure == 'List'):
            outEntries = self.lines

        # b) Entries organized into a dict (TO TEST LATER)
        elif (self.structure == 'Dict'):
            outEntries = [line for line in self.lines.values()]

        # c) Entries organized into a nested dict (TO IMPLEMENT LATER)
        #elif (self.structure == 'nestedDict'):

        ## Write entries into output bed file
        with open(outPath, 'w') as outFile:

            ## Write header
            fields = ['#ref', 'beg', 'end'] + optional

            if hasattr(outEntries[0], 'name'):
                fields.append('name')

            row = "\t".join(fields)
            outFile.write(row + '\n')

            ## Write entries 
            for entry in outEntries:

                fields = [entry.ref, str(entry.beg), str(entry.end)]

                # Add optional fields
                for field in optional:
                    fields.append(entry.optional[field])

                # Create row
                row = "\t".join(fields)
                outFile.write(row + '\n')

    def write_annovar(self, outPath):
        '''
        BED file writer. Write bed entries into a file format that can be used as input for annovar

        Input:
            1) outPath: path to output bed file
        '''
        ## Collect all the entries into a list
        # a) Entries organized into a list
        if (self.structure == 'List'):
            outEntries = self.lines

        # b) Entries organized into a dict (TO TEST LATER)
        elif (self.structure == 'Dict'):
            outEntries = [line for line in self.lines.values()]

        # c) Entries organized into a nested dict (TO IMPLEMENT LATER)
        #elif (self.structure == 'nestedDict'):

        ## Write entries into output bed file
        with open(outPath, 'w') as outFile:
            for entry in outEntries:

                fields = [entry.ref, str(entry.beg), str(entry.end), '0', '0']

                # Add name if available
                if 'name' in entry.optional:
                    fields.append('comments: ' + entry.optional['name'])

                # Create row
                row = "\t".join(fields)
                outFile.write(row + '\n')

    def organize_list(self, filePath, targetRefs):
        '''
        Organize bed file lines in a list

        Input:
            1) filePath: path to bed file
            2) targetRefs: list containing target references. Entries with refs not included in the list will not be loaded. 
                            If 'None' all the entries will be loaded
        
        Output:
            1) lines: list containing bed entries
        '''
        bedFile = open(filePath)
        lines = []
        header = []
        
        # For line in the file
        for line in bedFile:
            
            # Skip blank lines
            if not line:
                continue
            
            # Split data line into fields
            fields = line.split()                 

            # A) Header
            if line.startswith('#'):
                header = fields

            # B) Data line
            else:
                line = BED_entry(fields, header)

                if (targetRefs is None) or (line.ref in targetRefs):
                    lines.append(line)
        
        return lines

    def organize_dict(self, filePath, targetRefs):
        '''
        Organize bed file lines in a dictionary where each key will correspond to a reference and the corresponding value will be the list of lines in that reference
        
        Input:
            1) filePath: path to bed file
            2) targetRefs: list containing target references. Entries with refs not included in the list will not be loaded. 
                            If 'None' all the entries will be loaded

        Output:
            1) lines: dictionary containing bed entries
        '''

        bedFile = open(filePath)
        lines = {}
        header = []

        # For line in the file
        for line in bedFile:
            
            # Skip blank lines
            if not line:
                continue
            
            # Split data line into fields
            fields = line.split()       

            # A) Header
            if line.startswith('#'):
                header = fields

            # B) Data line
            else:
                line = BED_entry(fields, header)

                if (targetRefs is None) or (line.ref in targetRefs):

                    # a) Initialize reference and add line
                    if line.ref not in lines:
                        lines[line.ref] = [line]

                    # b) Add to preexisting reference
                    else:
                        lines[line.ref].append(line)

        return lines

    def organize_nestedDict(self, filePath, targetRefs):
        '''
        Organize bed file lines in a nested dictionary where first level keys correspond to the references and second level keys to bed entry names  

        Input:
            1) filePath: path to bed file
            2) targetRefs: list containing target references. Entries with refs not included in the list will not be loaded. 
                            If 'None' all the entries will be loaded

        Output:
            1) lines: nested dictionary containing bed entries

        NOTE: I need to update the code to take into account the bed header and optional fields
        '''

        bedFile = open(filePath)
        lines = {}
        header = []

        # For line in the file
        for line in bedFile:
            
            # Skip blank lines
            if not line:
                continue
            
            # Split data line into fields
            fields = line.split()       

            # A) Header
            if line.startswith('#'):
                header = fields

            # B) Data line
            else:
                line = BED_entry(fields, header)
                
                if (targetRefs is None) or (line.ref in targetRefs):

                    # A) Initialize reference and add line
                    if line.ref not in lines:
                        lines[line.ref] = {}
                        lines[line.ref][line.optional['name']] = [line]

                    # B) Add to preexisting reference
                    else:

                        # a) Initialize entry name
                        if line.optional['name'] not in lines[line.ref]:
                            lines[line.ref][line.optional['name']] = [line]

                        # b) Add to preexisting entry name
                        else:
                            lines[line.ref][line.optional['name']].append(line)

        return lines

    def group_entries_by_name(self):
        '''
        Organize bed file lines into a dictionary by feature name

        Output:
            1) groupedEntries: dictionary containing feature names as keys and as values a list of bed lines corresponding to each name
        '''
        groupedEntries = {}

        # For each bed entry
        for entry in self.lines:

            # a) Initialize entry name
            if entry.optional['name'] not in groupedEntries:
                groupedEntries[entry.optional['name']] = [entry]

            # b) Add to preexisting entry name
            else:
                groupedEntries[entry.optional['name']].append(entry)
                        
        return groupedEntries

class BED_entry():
    '''
    BED line class 
    '''
    number = 0 # Number of instances

    def __init__(self, fields, header):
        '''
        Initialize bed line

        Input:
            1. fields: list containing a bed feature (== data line)
            2. header: list containing bed header (required for parsing optional fields)
        '''
        BED_entry.number += 1 # Update instances counter
        self.id = 'BED_entry_' + str(BED_entry.number)

        ## Mandatory fields
        self.ref = str(fields[0])
        self.beg = int(fields[1])
        self.end = int(fields[2])
        self.clusterId = None

        ## Optional fields dictionary (Optional fields only considered if header provided)
        self.optional = {}

        for i in range(3, len(header), 1):
            self.optional[header[i]] = fields[i]       

class PAF():
    '''
    Class for dealing with files in PAF format. 
    '''
    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.alignments = []

    def read(self, filePath):
        '''
        PAF file reader. Read and store data line objects into a list:
        '''
        pafFile = open(filePath)

        # For line in the file
        for line in pafFile:
            
            # Skip comments and blank lines
            if line.startswith('#') or not line:
                continue

            fields = line.split() 

            ## Skip truncated lines
            if len(fields) < 12:
                print('[ERROR] Truncated line: ' + str(len(fields)) + ' fields')
                continue

            line = PAF_alignment(fields)
            self.alignments.append(line)

    def sortByLen(self):
        '''
        Sort alignments by query alignment length in descending order
        '''
        sortedAlignments = sorted(self.alignments, key=lambda x: x.alignmentLen(), reverse=True)
        return sortedAlignments

    ## [SR CHANGES]
    def sortNbMatches(self):
        '''
        Sort alignments by query alignment length
        '''
        sortedAlignments = sorted(self.alignments, key=lambda x: x.nbMatches, reverse=True)
        return sortedAlignments

    def chain(self, maxDist, maxPercOverlap):
        '''
        Chain PAF alignments based on alignment complementariety

        Input:
            1. maxDist: maximum distance between both ranges 
            2. maxPercOverlap: maximum percentage of overlap between ranges

        Output:
            1. chain: PAF_chain object instance
        '''
        ## 1. Sort alignments by decreasing query alignment length 
        sortedAlignments = self.sortByLen()

        ## 2. Pick longest alignment and initiate chain
        longest = sortedAlignments[0]
        chain = PAF_chain([longest])
        
        # remove alignment 
        del sortedAlignments[0]

        roundCounter = 1

        ## 3. Attemp to extend the chain with complementary alignments
        while True:

            # START ALIGNMENT CHAIN EXTENSION ROUND
            # Initialize boolean as not complementary alignment found
            complBool = False

            # Go through all the available alignments 
            for index, alignment in enumerate(sortedAlignments):

                ## Assess if alignment complementary to the chain
                chainBeg, chainEnd = chain.interval()
                complBool, orientation = gRanges.complementary(chainBeg, chainEnd, alignment.qBeg, alignment.qEnd, maxDist, maxPercOverlap)

                ## Complementary alignment found 
                if complBool:

                    ## Add alignment to the chain 
                    # a) Add to the chain begin
                    if orientation == "LEFT":
                        chain.alignments.insert(0,alignment)

                    # b) Add to the chain end
                    else:
                        chain.alignments.append(alignment)

                    ## Remove from list
                    del sortedAlignments[index]

                    ## Stop once complementary found
                    break

            roundCounter += 1

            # STOP CHAIN EXTENSION IF: 
            # a) No complementary alignment found in the last round OR
            # b) Mo alignments left
            if complBool == False or not sortedAlignments:
                break

        return chain

    def hits2dict(self):
        '''
        Reorganize hits into a dictionary
        '''

        hitsDict = {}

        # For each hit
        for hit in self.alignments:

            # Initialize list
            if hit.qName not in hitsDict:
                hitsDict[hit.qName] = []
            
            # Add hit to list
            hitsDict[hit.qName].append(hit)
    
        return hitsDict

class PAF_alignment():
    '''
    PAF entry class 
    '''
    number = 0 # Number of instances

    def __init__(self, fields):
        '''
        Initialize paf line
        '''
        PAF_alignment.number += 1 # Update instances counter
        self.id = 'PAF_alignment_' + str(PAF_alignment.number)
        self.qName = str(fields[0])
        self.qLen = int(fields[1])
        self.qBeg = int(fields[2])
        self.qEnd = int(fields[3])
        self.strand = str(fields[4])
        self.tName = str(fields[5])
        self.tLen = int(fields[6])
        self.tBeg = int(fields[7])
        self.tEnd = int(fields[8])
        self.nbMatches = int(fields[9])
        self.blockLen = int(fields[10])
        self.MAPQ = int(fields[11])   

    def alignmentString(self):
        '''
        Create string with alignment information: tBeg-tEnd,alignmentLen(based on the query)
        '''
        string = str(self.alignmentLen()) + ',' + str(self.tBeg) + '-' + str(self.tEnd)
        return string
    
    def alignmentLen(self):
        '''
        Compute the query alignment length
        '''

        return self.qEnd - self.qBeg

    def alignmentPerc(self):
        '''
        Compute the query alignment length percentage
        '''
        percLen = float(self.alignmentLen()) / self.qLen * 100

        return percLen

class PAF_chain():
    '''    
    Chain of complementary PAF alignments  
    '''
    def __init__(self, alignments):
        '''
        Initialize chain instance. 
        
        Input:
            1. alignments. List of PAF_alignment instances
        '''
        self.alignments = alignments

    def interval(self):
        '''
        Return query interval covered by the chain
        '''
        firstAlignment = self.alignments[0]
        lastAlignment = self.alignments[-1]
        
        return firstAlignment.qBeg, lastAlignment.qEnd

    def interval_template(self):
        '''
        Return query interval covered by the chain
        '''
        firstAlignment = self.alignments[0]
        lastAlignment = self.alignments[-1]
        
        return firstAlignment.tBeg, lastAlignment.tEnd

    def perc_query_covered(self):
        '''
        Compute the percentage of the query sequence covered by the chain of alignments
        '''
        # a) No alignments available
        if len(self.alignments) == 0:
            percCovered = 0

        # b) Alignments available
        else:

            ## Compute the number of bases covered
            beg, end = self.interval()
            alignmentLen = end - beg

            ## Compute the percentage of bases covered
            percCovered = float(alignmentLen)/self.alignments[0].qLen*100

        return percCovered

    def nb_templates(self):
        '''
        Compute the number of different templates in the chain
        '''
        templates = list(set([alignment.tName for alignment in self.alignments]))
        nbTemplates = len(templates)

        return nbTemplates, templates

class VCF():
    '''
    VCF class
    '''
    def __init__(self):
        '''
        '''
        self.header = None
        self.variants = []  # List of variants
        self.info_order = []
        self.format_order = []

    def read(self, filePath):
        '''
        Read VCF file
        '''
        vcfFile = open(filePath)

        ## 1. Read VCF and differenciate between header and variant entries
        header = []
        variants = []

        # For line in the file
        for line in vcfFile:
            
            # Skip blank lines
            if not line:
                continue    

            # Remove trailing spaces and new line
            line = line.rstrip()

            # a) Header
            if line.startswith('#'):
                header.append(line)

            # b) Variant
            elif not line.startswith('#'):
                variants.append(line)
        
        ## 2. Read header entries
        self.read_header(header)

        ## 3. Read variant entries
        self.read_variants(variants)

    def read_header(self, header):
        '''
        Read VCF file header

        Input:
            1. header: list of header lines
        '''
        species = ''
        refLengths = {} 
        info = {}
        gtFormat = {}

        source = None
        build = None
        sampleIds = None

        ## Read header line by line
        for line in header:

            ## A) Source
            if line.startswith('##source'):
                source = line.split('=', 1)[1]
            
            ## B) Reference
            elif line.startswith('##reference'):
                build = line.split('=', 1)[1]

            ## C) Contig
            elif line.startswith('##contig'):
                substring = re.search('<(.*)>', line)

                for field in substring.group(1).split(','):
                    key, value = field.split('=', 1)

                    if key == 'ID':
                        cId = value

                    if key == 'length':
                        cLen = value

                    if key == 'species':
                        species = value

                refLengths[cId] = cLen
            
            ## D) Info
            elif line.startswith('##INFO'):
                substring = re.search('<(.*)>', line)
                values = []

                for field in substring.group(1).split(','):
                    fields = field.split('=', 1)
                    
                    if len(fields) == 2:
                        values.append(fields[1])

                info[values[0]] = values[1:]
                self.info_order.append(values[0])

            ## E) Genotype format
            elif line.startswith('##FORMAT'):
                substring = re.search('<(.*)>', line)
                values = []

                for field in substring.group(1).split(','):
                    fields = field.split('=', 1)
                    
                    if len(fields) == 2:
                        values.append(fields[1])

                gtFormat[values[0]] = values[1:]
                self.format_order.append(values[0])

            ## F) Colnames line
            elif line.startswith('#CHROM'):

                fields = line.split("\t")

                ## a) Multisample VCF file
                if len(fields) > 8:
                    sampleIds = [fields[i] for i in range(9, len(fields))]
                
        ## Create VCF header object
        self.create_header(source, build, species, refLengths, info, gtFormat, sampleIds)

    def read_variants(self, entries):
        '''
        Read variant entries from VCF file

        Input:
            1. entries: list VCF variant entries
        '''
        # For each entry
        for entry in entries:

            fields = entry.split("\t")

            ## Parse info field
            INFO = {}
            infoFields = fields[7].split(';')

            # For each feature at info
            for feature in infoFields:
                feature = feature.split('=')

                # a) Flag 
                if len(feature) == 1:
                    key = str(feature[0])
                    value = True
                    INFO[key] = value

                # b) Key and value pair
                else:
                    key = feature[0]
                    value = feature[1]
                    INFO[key] = value

            ## Parse sample genotypes if multisample VCF
            FORMAT = {}

            if len(fields) > 8:
                
                samplesGt = [fields[i] for i in range(9, len(fields))]   

                ## For each sample Genotype filed
                for sampleGt, sampleId in zip(samplesGt, self.header.sampleIds):
                    FORMAT[sampleId] = {}
                    gtFields = sampleGt.split(':')
                    
                    for key, value in zip(self.format_order, gtFields):
                        FORMAT[sampleId][key] = value

            ## Create VCF variant object
            fields = fields[0:7] + [INFO] + [FORMAT]
            variant = VCF_variant(fields)

            ## Add variant to the VCF
            self.add(variant)

    def add(self, variant):
        '''
        Add VCF_variant instance to the VCF 
        '''
        self.variants.append(variant)

    def create_header(self, source, build, species, refLengths, info, gtFormat, sampleIds):
        '''
        Create VCF header

        Input:
            1. source: Software version used to generate the insertion calls
            2. build: Reference genome build
            3. species: Specie
            4. refLengths: Dictionary containing reference ids as keys and as values the length for each reference
            5. info: Dictionary containing data to include at INFO. Dictionary keys will correspond
                     to INFO entry identifiers while values will be 3 element lists corresponding to Number, Type and 
                     Description fields for an INFO entry.
            6. gtFormat: Dictionary containing data to include at FORMAT field. Dictionary keys will correspond 
                     to FORMAT entry identifiers while values will be 3 element lists corresponding to Number, Type and 
                     Description fields for an FORMAT entry.
            7. sampleIds: list of sample identifiers. None if not multisample VCF file

        Output:
            Create and include header object at VCF class
        '''
        self.header = VCF_header(source, build, species, refLengths, info, gtFormat, sampleIds)

    def sort(self):
        '''
        Sort variants in increasing coordinates ordering
        '''
        self.variants = sorted(self.variants, key=lambda variant: (variant.chrom, variant.pos))

    def write(self, infoIds, formatIds, outName, outDir):
        '''
        Write VCF into output file

        Input:
            1. infoIds: Array of info fields to be listed (same order as the list)
            2. formatIds: Array of format fields to be listed (same order as the list)
            2. outName: Output file name
            3. outDir: Output directory

        Output: Write VCF file 
        '''     
        ## 1. Open output filehandle    
        outFile = outDir + '/' + outName + '.vcf'
        outFile = open(outFile, 'w')

        ## 2. Write header
        header = self.header.build_header(infoIds, formatIds)
        outFile.write(header)

        ## 3. Write variants
        for variant in self.variants:

            INFO = variant.build_info(infoIds)

            # a) Regular VCF
            if not formatIds:
                row = "\t".join([variant.chrom, str(variant.pos), str(variant.ID), variant.ref, variant.alt, variant.qual, variant.filter, INFO, "\n"])

            # b) Multi-sample VCF
            else:
                genotypes = variant.build_genotypes(formatIds, self.header.sampleIds)
                row = "\t".join([variant.chrom, str(variant.pos), str(variant.ID), variant.ref, variant.alt, variant.qual, variant.filter, INFO, ':'.join(formatIds)] + genotypes + ["\n"])
            
            outFile.write(row)

        ## Close output file
        outFile.close()

    def ins2fasta(self, itype, fam, outDir):
        '''
        Write inserted sequence into a fasta

        Input:
            1. itype: Target insertion type
            2. fam: Target family
            3. outDir: Output directory

        Output: 
            1. fasta: path to fasta file containing inserted sequences
        '''     
        ## 1. Initialize fasta object
        fasta = FASTA() 

        ## 2. Collect inserted sequences
        for variant in self.variants:

            if (variant.info['ITYPE'] == itype) and ('FAM' in variant.info) and (variant.info['FAM'] == fam):
                
                seqId = variant.chrom + '_' + str(variant.pos) + '_' + variant.info['FAM'] 
                seq = variant.info['INSEQ'] 

                fasta.seqDict[seqId] = seq

        ## 3. Write inserted sequences into fasta file
        outFile = outDir + '/' + itype + '_' + fam + '.fa'
        fasta.write(outFile)

        return outFile

class VCF_header():
    '''
    VCF header class
    '''
    
    def __init__(self, source, build, species, refLengths, info, gtFormat, sampleIds):
        '''
        Initialize VCF header
        '''
        self.source = source
        self.build = build
        self.species = species
        self.refLengths = refLengths
        self.info = info
        self.format = gtFormat
        self.sampleIds = sampleIds

    def build_header(self, infoIds, formatIds):
        '''
        Build VCF header string
        
        Input:
            1. infoIds: Array of info fields to be listed (same order as the list)
            2. formatIds: Array of format fields to be listed (same order as the list)

        Output:
            1. header: Header string
        '''
        ## 1. Collect general features
        date = time.strftime("%Y%m%d")

        data = {
            'date': date,
            'source': self.source,
            'reference': self.build
        }

        template = """##fileformat=VCFv4.2\n##fileDate={date}\n##source={source}\n##reference={reference}\n"""
        general = template.format(**data)
        
        ## 2. Build contigs
        contigs = self.build_contigs()

        ## 3. Build info
        info = self.build_info(infoIds)
        
        ## 3. Build format
        formatGt = self.build_format(formatIds)

        ## 5. Column data names
        # a) Regular VCF
        if not self.sampleIds:
            colnames = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', "\n"])

        # b) Multi-sample VCF
        else:
            colnames = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + self.sampleIds + ["\n"])


        ## 6. Join all the info
        header = ''.join([general, contigs, info, formatGt, colnames])
        return header

    def build_contigs(self):
        '''
        Build header INFO string 

        Output:
            1. info: Header CONTIG string
        '''

        entries = []
        for refId in sorted(self.refLengths.keys()):

            data = {
                'refId': refId,
                'build': self.build,
                'length': self.refLengths[refId],
                'species': self.species,
            }
            
            template = """##contig=<ID={refId},assembly={build},length={length},species={species}>\n"""
            entry = template.format(**data)
            entries.append(entry)

        CONTIGS = ''.join(entries)
        
        return CONTIGS

    def build_info(self, IDS):
        '''
        Build header INFO string 

        Input:
            1. IDS: Array of info fields to be listed (same order as the list)

        Output:
            1. INFO: Header INFO string
        '''
        entries = []

        for ID in IDS:

            data = {
                'ID': ID,
                'Number': self.info[ID][0],
                'Type': self.info[ID][1],
                'Description': self.info[ID][2]
            }
            template = """##INFO=<ID={ID},Number={Number},Type={Type},Description=\"{Description}\">\n"""
            entry = template.format(**data)
            entries.append(entry)

        INFO = ''.join(entries)

        return INFO

    def build_format(self, IDS):
        '''
        Build header FORMAT string 

        Input:
            1. IDS: Array of FORMAT fields to be listed (same order as the list)

        Output:
            1. FORMAT: Header FORMAT string
        '''
        entries = []

        for ID in IDS:

            data = {
                'ID': ID,
                'Number': self.format[ID][0],
                'Type': self.format[ID][1],
                'Description': self.format[ID][2],
            }
            
            template = """##FORMAT=<ID={ID},Number={Number},Type={Type},Description={Description}>\n"""
            entry = template.format(**data)
            entries.append(entry)

        FORMAT = ''.join(entries)

        return FORMAT

class VCF_variant():
    '''
    VCF variant class
    '''
    number = 0 # Number of instances

    def __init__(self, fields):
        '''
        Initialize VCF variant class
        '''
        VCF_variant.number += 1 # Update instances counter
        self.id = 'variant_' + str(VCF_variant.number)

        self.chrom = fields[0]
        self.pos = int(fields[1])
        self.ID = fields[2]
        self.ref = fields[3]
        self.alt = fields[4]
        self.qual = fields[5]
        self.filter = fields[6]
        self.info = fields[7]
        self.format = fields[8]
        self.clusterId = None

    def build_info(self, IDS):
        '''
        Create info field string
        
        Input:
            1. IDS: Array of info fields to be listed (same order as the list)

        Output:
            1. INFO: Info string
        '''        
        INFO = ''

        for index, ID in enumerate(IDS):
            
            ## Include field
            if (ID in self.info) and (self.info[ID] is not None) and (self.info[ID] is not False):

                # a) Boolean
                if isinstance(self.info[ID], bool):
                    entry = ID

                # b) Not boolean
                else:
                    entry = ID + '=' + str(self.info[ID])

                INFO = INFO + entry + ';'

        INFO = INFO[:-1]
        return INFO

    def build_genotypes(self, IDS, sampleIds):
        '''
        Create format field string
        
        Input:
            1. IDS: Array of format fields to be listed (same order as the list)

        Output:
            1. genotypes: list
        '''      
        genotypes = []

        for sampleId in sampleIds:
                    
            genotype = ':'.join([self.format[sampleId][ID] for ID in IDS if ID in self.format[sampleId]])
            genotypes.append(genotype)

        return genotypes

    def pos_interval(self):
        '''
        Compute position interval
        '''        

        if 'CIPOS' in self.info:
            ciBeg, ciEnd = self.info['CIPOS'].split(',')
            beg = self.pos + int(ciBeg)
            end = self.pos + int(ciEnd)
        else:
            beg = self.pos
            end = self.pos

        return beg, end

class INS_cluster():
    '''
    '''
    number = 0 # Number of instances

    def __init__(self, variants):
        '''
        Initialize class instance
        '''
        ## Set cluster id
        INS_cluster.number += 1 # Update instances counter
        self.id = 'CLUSTER_' + str(INS_cluster.number)

        # Define list of events composing the cluster 
        self.variants = variants

        # Update event's clusterId attribute
        for event in self.variants:
            event.clusterId = self.id   

    def add(self, variants2add):
        '''
        Incorporate variants into the cluster 

        Input:
            1. variants2add: List of VCF variant objects to be added to the cluster
        '''
        # Add events to the cluster  
        self.variants = self.variants + variants2add

    def consensus(self):
        '''
        Select representative INS variant for the cluster. If there are hifi variants, 
        choose one arbitraty HiFi as representative

        Note: right now I will pick an arbitrary variant. Improve in the future. Ideas
            - Improve consensus sequence quality by making a consensus of consensus sequences
        '''

        ## Select arbitrary event as consensus (improve later)
        platforms = sorted(set([variant.platform for variant in self.variants]))

        # a) Variant identified in HiFi 
        if 'HIFI' in platforms:
            consensus = [variant for variant in self.variants if variant.platform == 'HIFI'][0]
                
        # b) Variant not identified in HiFi 
        else:
            consensus = self.variants[0]

        consensus.info['LEAD_SAMPLE'] = consensus.sampleId

        ## Create platform info field
        consensus.info['PLATFORMS'] = ','.join(platforms)

        ## Add list of samples where the event is identified
        sampleIds = ','.join(set([event.sampleId for event in self.variants]))
        consensus.info['SAMPLES'] = sampleIds

        return consensus

class PSL():
    '''
    Class for dealing with files in PSL format. 
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.alignments = []

    def read(self, filePath):
        '''
        PSL file reader. Read and store data line objects into a list:
        '''
        fh = open(filePath)

        # For line in the file
        for line in fh:
            
            # Skip comments and blank lines
            if line.startswith('#') or not line:
                continue

            fields = line.split() 

            self.alignments.append(PSL_alignment(fields))

    def alignments2dict(self):
        '''
        Organize the alignments into a dictionary by query name
        '''
        alignDict = {}

        for alignment in self.alignments:

            if alignment.qName not in alignDict:
                alignDict[alignment.qName] = []

            alignDict[alignment.qName].append(alignment)

        return alignDict

    def hits2clipping(self, clippings):
        '''
        Add hits to clipping events as supplementary alignments.
        '''
        # For each alignment
        for alignment in self.alignments:

            ## Compute CIGAR
            clipping = clippings[alignment.qName]
            hitLen = alignment.qEnd - alignment.qBeg
            hardClipLen = len(clipping.readSeq) - hitLen

            if clipping.clippedSide == 'left':
                CIGAR = str(hitLen) + 'M' + str(hardClipLen) + 'H'
            
            else:
                CIGAR = str(hardClipLen) + 'H' + str(hitLen) + 'M' 

            ## Create SA string
            fields = [alignment.tName, str(alignment.tBeg), alignment.strand, CIGAR, '60', str(alignment.misMatches)] 
            SA = ','.join(fields) + ';'

            ## Add SA to the clipping event
            if clippings[alignment.qName].supplAlignment is None:
                clippings[alignment.qName].supplAlignment = SA
            else:
                clippings[alignment.qName].supplAlignment = clippings[alignment.qName].supplAlignment + SA 

    def filter_align_perc(self, minPerc):
        '''
        Filter the alignments by the percentage of their sequence that is aligned

        Input:
            1. minPerc: minimum percentage of aligned sequence
        
        Output:
            1. filteredAlignments: list of filtered alignments
        '''
        filteredAlignments = []

        ## For each alignment
        for alignment in self.alignments:

            ## Select alignment as enough perc of sequence aligned
            if alignment.perc_query_covered() >= minPerc:

                filteredAlignments.append(alignment)

        return filteredAlignments

    def filter_nb_hits(self, maxNbHits):
        '''
        Filter alignments based on their number of hits. 

        Input:
            1. maxNbHits: maximum number of hits allowed

        Output.
            1. filteredAlignments: list of filtered alignments
        '''
        filteredAlignments = []

        ## 1. Organize alignments by query name
        alignDict = self.alignments2dict()

        ## 2. Filter out alignments for queries with more than X hits
        filteredAlignments = []

        for queryName in alignDict:
            nbHits = len(alignDict[queryName])

            if nbHits <= maxNbHits:
                filteredAlignments = filteredAlignments + alignDict[queryName]

        return filteredAlignments

class PSL_alignment():
    '''
    PSL alignment class
    '''

    def __init__(self, fields):
        '''
        '''

        # Define blat alignment variables
        self.matches = int(fields[0])
        self.misMatches = int(fields[1])
        self.repMatches = int(fields[2])
        self.nCount = int(fields[3])
        self.qNumInsert = int(fields[4])
        self.qBaseInsert = int(fields[5])
        self.tNumInsert = int(fields[6])
        self.tBaseInsert = int(fields[7])
        self.strand = fields[8]
        self.qName = fields[9]
        self.qSize = int(fields[10])
        self.qBeg = int(fields[11])
        self.qEnd = int(fields[12])
        self.tName = fields[13]
        self.tSize = int(fields[14])
        self.tBeg = int(fields[15])
        self.tEnd = int(fields[16])
        self.blockCount = int(fields[17])
        self.blockSizes = fields[18]
        self.qStarts = fields[19]
        self.tStarts = fields[20]

        # Other
        self.alignType = ""

    def perc_query_covered(self):
        '''
        Compute percentage of query aligned 
        '''
        percQueryCovered = float(self.qEnd - self.qBeg) / self.qSize * 100  
        return percQueryCovered

    def rev_complement(self):
        '''
        Make the reverse complementary aligment. This would be the alignment information of the reverse complementary original sequence
        '''
        ## Reverse complement strand
        switchStrandDict = {'+': '-', '-': '+'}
        self.strand = switchStrandDict[self.strand]

        ## Reverse complement query start and end positions
        updatedBeg = self.qSize - self.qEnd
        updatedEnd = self.qSize - self.qBeg

        self.qBeg = updatedBeg
        self.qEnd = updatedEnd

def pslQueryRefDict(pslPath):
    '''
    Read BLAT results and store qName and tName in a dictionary

    Input:
        1. pslPath: path to blat result file (psl format)
    Output:
        1. pslDict: dictionary -> pslDict[qName] = tName 
    '''
    # Read PSL
    pslClipping = PSL()
    pslClipping.read(pslPath)
    ## TODO SR: mirar filtros
    pslDict = {}
    for pslAlign in pslClipping.alignments:
        if pslAlign.qName in pslDict.keys():
            pslDict[pslAlign.qName].append(pslAlign.tName)
        else:
            pslDict[pslAlign.qName] = []
            pslDict[pslAlign.qName].append(pslAlign.tName)
    return pslDict


class mosaiCatcherSVs():
    '''
    Class for dealing with mosaiCatcher output file containing structural variants (SV)
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.SVs = None
        self.structure =  None

    def read(self, filePath, structure, targetRefs):
        '''
        File reader. Read and store SV calls into a data structure

        Input:
            1) filePath: path to MosaiCatcher output file
            2) structure: data structure where SV calls will be stored. 3 Possible structures:
                - 'List': SVs organized into a list
                - 'Dict': SVs organized into a dictionary where each key will correspond to a reference and the corresponding value will be the list of SVs in that reference
                - 'NestedDict': SVs organized into a nested dictionary with first level keys corresponding to the cell, second level to the reference and the corresponding value the list of SVs
            3) targetRefs: list containing target references. SVs with refs not included in the list will not be loaded. 
                            If 'None' all the SVs will be loaded

        Initialize SVs attribute as output
        '''
        self.structure = structure

        # a) Organize SVs into a list
        if (self.structure == 'List'):

            self.SVs = self.organize_list(filePath, targetRefs)

        # b) Organize BED entries into a dict
        elif (self.structure == 'Dict'):
            
            self.SVs = self.organize_dict(filePath, targetRefs)

        # b) Organize BED entries into a nested dict
        elif (self.structure == 'NestedDict'):
            
            self.SVs = self.organize_nestedDict(filePath, targetRefs)

        # d) Unkown data type structure provided
        else:
            print('[ERROR] MosaiCatcher SV file reader. Unknown structure provided: ', structure)
            sys.exit(1)

    def organize_list(self, filePath, targetRefs):
        '''
        Organize bed file lines in a list

        Input:
            1) filePath: path to MosaiCatcher output file
            2) targetRefs: list containing target references. SVs with refs not included in the list will not be loaded. 
                            If 'None' all the SVs will be loaded
        
        Output:
            1) SVs: list containing SVs
        '''
        mosaiCatcherOut = open(filePath)
        SVs = []
        
        # For line in the file
        for line in mosaiCatcherOut:
            
            # Skip blank lines and header
            if (not line) or line.startswith('chrom'):
                continue
            
            # Split data line into fields
            fields = line.split()   

            # Create SV object              
            SV = SV_call(fields)

            # Add SV to the list
            if (targetRefs is None) or (SV.ref in targetRefs):
                SVs.append(SV)
        
        return SVs

    def organize_dict(self, filePath, targetRefs):
        '''
        Organize bed file lines in a dictionary where each key will correspond to a reference and the corresponding value will be the list of lines in that reference
        
        Input:
            1) filePath: path to MosaiCatcher output file
            2) targetRefs: list containing target references. SVs with refs not included in the list will not be loaded. 
                            If 'None' all the SVs will be loaded

        Output:
            1) SVs: dictionary containing SVs. One key per reference and the list of SVs per reference as value
        '''
        mosaiCatcherOut = open(filePath)
        SVs = {}

        # For line in the file
        for line in mosaiCatcherOut:
            
            # Skip blank lines and header
            if (not line) or line.startswith('chrom'):
                continue
            
            # Split data line into fields
            fields = line.split()       
            SV = SV_call(fields)

            if (targetRefs is None) or (SV.ref in targetRefs):

                # a) Initialize reference and add SV
                if SV.ref not in SVs:
                    SVs[SV.ref] = [SV]

                # b) Add to preexisting reference
                else:
                    SVs[SV.ref].append(SV)

        return SVs

    def organize_nestedDict(self, filePath, targetRefs):
        '''
        Organize SVs into a nested dictionary with first level keys corresponding to the cell, second level to the reference and the corresponding value the list of SVs
        
        Input:
            1) filePath: path to MosaiCatcher output file
            2) targetRefs: list containing target references. SVs with refs not included in the list will not be loaded. 
                           If 'None' all the SVs will be loaded

        Output:
            1) SVs: dictionary containing SVs. One key per cellId and reference and the lists of SVs as values
        '''
        mosaiCatcherOut = open(filePath)
        SVs = {}

        # For line in the file
        for line in mosaiCatcherOut:
            
            # Skip blank lines and header
            if (not line) or line.startswith('chrom'):
                continue
            
            # Split data line into fields
            fields = line.split()       
            SV = SV_call(fields)

            if (targetRefs is None) or (SV.ref in targetRefs):

                # A) Initialize cell 
                if SV.cellId not in SVs:
                    SVs[SV.cellId] = {}
                    SVs[SV.cellId][SV.ref] = [SV]

                # B) Add to preexisting cell
                else:

                    # a) Initialize reference and add SV
                    if SV.ref not in SVs[SV.cellId]:
                        SVs[SV.cellId][SV.ref] = [SV]

                    # b) Add to preexisting reference
                    else:
                        SVs[SV.cellId][SV.ref].append(SV)
                        
        return SVs

class SV_call():
    '''
    SV call class 
    '''
    number = 0 # Number of instances

    def __init__(self, fields):
        '''
        Initialize paf line
        '''
        SV_call.number += 1 # Update instances counter
        self.id = 'SV_call_' + str(SV_call.number)

        ## Coordinates
        self.ref = str(fields[0])
        self.ref = self.ref.replace('chr', '') # Remove chr string from the identifier
        self.beg = int(float(fields[1]))
        self.end = int(float(fields[2]))
        self.dist2end = None                   # Distance to the telomere end
        
        ## Identifiers
        self.sampleId = str(fields[3])
        self.cellId = str(fields[4])

        ## Strand state
        self.state = str(fields[5])  # WW, CC, WC  
        self.scalar = int(fields[6]) # ?? seems to be always one
        self.nbBins = int(fields[7]) # Number of genomic bins spanned by the SV call

        ## Structural variation (SV) type
        str(fields[8]).split('_')

        self.SV_type1, self.haplo1 = self.parse_sv_call_name(fields[8])  # Most likely SV and haplotype   
        self.haploCode1 = str(fields[9]) # Haplotype code for 1. Unique identifier for SV type and haplotype

        self.SV_type2, self.haplo2 = self.parse_sv_call_name(fields[10]) # Second most likely SV and haplotype   
        self.haploCode2 = str(fields[11]) # Haplotype code for 2. Unique identifier for SV type and haplotype

        ## Statistics
        self.llr2ref = float(fields[12]) # Log likelyhood ratio between most likely SV and the reference
        self.llr2SV_type2 = float(fields[13]) # Log likelyhood ratio between most likely SV and second
        self.AF = float(fields[14]) # Variant allele fraction (in the population of cells)

    def len(self):
        '''
        Compute SV length

        Output:
            1. length: structural variant length
        '''
        length = self.end - self.beg

        return length
    
    def parse_sv_call_name(self, string):
        '''
        Parse SV call name string

        Input:
            1. string: 'sv_call_name' or 'sv_call_name_2nd' fields from MosaiCatcher SVs output file

        Output:
            1. SV_type: Structural variant type (del, dup, idup, inv, complex)
            2. haplotype: Haplotype (h1 or h2)
        '''
        fields = string.split("_")   

        # a) Haplotype info not available
        if len(fields) == 1:
            SV_type = fields[0]
            haplotype = None           

        # b) Haplotype info available
        else:
            SV_type = fields[0]
            haplotype = fields[1]           

        return SV_type, haplotype

class TRF():
    '''
    Class for dealing with Tandem Repeat Finder output files
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.callsDict = {}

    def read(self, filePath):
        '''
        TRF file reader. Read and store tandem repeat calls a dictionary:
        '''
        trfFile = open(filePath)

        # For line in the file
        for line in trfFile:
            
            ## Skip comments and blank lines
            if line.startswith('#') or not line:
                continue

            ## Read sequence id
            if line.startswith('@'):
                seqId = line[1:-1]

                ## Initialize object
                if seqId not in self.callsDict:
                    self.callsDict[seqId] = TRF_calls(seqId)
            
                continue

            ## Create tandem repeat object and add to the objects dictionary
            self.callsDict[seqId].add(line)

    def add_seq_size(self, FASTA):
        '''
        Add the size to every analyzed sequence for tandem repeats
        '''
        for seqId in FASTA.seqDict:

            seqLen = len(FASTA.seqDict[seqId])

            if seqId in self.callsDict:
                self.callsDict[seqId].seqLen = seqLen

                for call in self.callsDict[seqId].calls:
                    call.seqLen = seqLen


class TRF_calls():
    '''
    TRF calls for a given sequence
    '''

    def __init__(self, seqId):
        '''
        Initialize empty class instance
        '''
        self.seqId = seqId
        self.seqLen = None
        self.calls = []

    def add(self, line):
        '''
        Add tandem repeat call
        '''
        ## Create TRF call and add to the list
        self.calls.append(TRF_call(line)) 

    def sortByLen(self):
        '''
        Sort alignments by query alignment length in descending order
        '''
        sortedLen = sorted(self.calls, key=lambda x: x.trLen, reverse=True)
        return sortedLen

    def is_tandem_repeat(self, minPerc, maxDist, maxPercOverlap):
        '''
        Assess if the sequence is a tandem repeat or not based on the % explained by TRF calls

        NOTE: Now I will use a very simple criteria for starting that will rely on picking the longest call and discarding the others. Later look into
        more sophisticated approaches which will search for complementary TR calls
        '''
        ## 0. Abort assessment if seqlen not available 
        if self.seqLen is None:
            return None, None, None, None

        ## 1. Chain TRF hits
        chain = self.chain(maxDist, maxPercOverlap)

        ## 2. Percentage sequence resolved based on the call
        fraction = chain.perc_query_covered()

        ## 3. Make decision
        if fraction >= minPerc:
            is_tandem_repeat = True

            ## Collect info
            # - Number of Motifs
            # - Motif sequences
            nbMotifs = len(chain.alignments)
            motifs = ','.join([hit.motifSeq for hit in chain.alignments])

        else:
            is_tandem_repeat = False
            nbMotifs = None
            motifs = None


        return is_tandem_repeat, fraction, nbMotifs, motifs

    def SVA_hexamer(self, seq, strand, maxDist, maxPercOverlap):
        '''
        Assess if the sequence is a tandem repeat or not based on the % explained by TRF calls

        NOTE: Now I will use a very simple criteria for starting that will rely on picking the longest call and discarding the others. Later look into
        more sophisticated approaches which will search for complementary TR calls
        '''
        ### 1. Select hexamer hits (CT rich repeats)
        hits = []

        for hit in self.calls:

            dist2beg = hit.beg
            dist2end = len(seq) - hit.end

            if (strand == '+') and (dist2beg < 300) and (hit.percC > 50) and (hit.percT > 20) and (len(hit.trSeq) > 10):
                hits.append(hit)
            
            elif (strand == '-') and (dist2end < 300) and (hit.percG > 50) and (hit.percA > 20) and (len(hit.trSeq) > 10):
                hits.append(hit)

        ### 2. Make hexamer call
        # a) No hit
        if not hits:
            hexamerSeq = None
            trimmedSeq = seq

        # b) Candidate hits
        else:
            ## Replace all hits by only hexamere hits
            self.calls = hits

            ## Chain hexamer hits
            chain = self.chain(maxDist, maxPercOverlap)
            hexBeg, hexEnd = chain.interval()
            hexamerSeq = seq[hexBeg:hexEnd]
            trimmedSeq = seq[:hexBeg] + seq[hexEnd:]

        return hexamerSeq, trimmedSeq

    def chain(self, maxDist, maxPercOverlap):
        '''
        Chain PAF alignments based on alignment complementariety

        Input:
            1. maxDist: maximum distance between both ranges 
            2. maxPercOverlap: maximum percentage of overlap between ranges

        Output:
            1. chain: PAF_chain object instance
        '''
        ## 1. Sort alignments by decreasing query alignment length 
        sortedAlignments = self.sortByLen()

        ## 2. Pick longest alignment and initiate chain
        longest = sortedAlignments[0]
        chain = TRF_chain([longest])
        
        # remove alignment 
        del sortedAlignments[0]

        roundCounter = 1

        ## 3. Attemp to extend the chain with complementary alignments
        while True:

            # START ALIGNMENT CHAIN EXTENSION ROUND
            # Initialize boolean as not complementary alignment found
            complBool = False

            # Go through all the available alignments 
            for index, alignment in enumerate(sortedAlignments):

                ## Assess if alignment complementary to the chain
                chainBeg, chainEnd = chain.interval()
                complBool, orientation = gRanges.complementary(chainBeg, chainEnd, alignment.beg, alignment.end, maxDist, maxPercOverlap)

                ## Complementary alignment found 
                if complBool:

                    ## Add alignment to the chain 
                    # a) Add to the chain begin
                    if orientation == "LEFT":
                        chain.alignments.insert(0,alignment)

                    # b) Add to the chain end
                    else:
                        chain.alignments.append(alignment)

                    ## Remove from list
                    del sortedAlignments[index]

                    ## Stop once complementary found
                    break

            roundCounter += 1

            # STOP CHAIN EXTENSION IF: 
            # a) No complementary alignment found in the last round OR
            # b) Mo alignments left
            if complBool == False or not sortedAlignments:
                break

        return chain

class TRF_chain():
    '''    
    Chain of complementary PAF alignments  
    '''
    def __init__(self, alignments):
        '''
        Initialize chain instance. 
        
        Input:
            1. alignments. List of PAF_alignment instances
        '''
        self.alignments = alignments

    def interval(self):
        '''
        Return query interval covered by the chain
        '''
        firstAlignment = self.alignments[0]
        lastAlignment = self.alignments[-1]
        
        return firstAlignment.beg, lastAlignment.end

    def perc_query_covered(self):
        '''
        Compute the percentage of the query sequence covered by the chain of alignments
        '''
        # a) No alignments available
        if len(self.alignments) == 0:
            percCovered = 0

        # b) Alignments available
        else:

            ## Compute the number of bases covered
            beg, end = self.interval()
            alignmentLen = end - beg

            ## Compute the percentage of bases covered
            percCovered = float(alignmentLen)/self.alignments[0].seqLen*100

        return percCovered

class TRF_call():
    '''
    TRF call for a given sequence
    '''
    def __init__(self, line):
        '''
        Initialize empty class instance
        '''
        # Split data line into fields
        fields = line.split()        

        # Create attributes
        self.beg = int(fields[0]) - 1
        self.end = int(fields[1]) - 1
        self.repeat_size = float(fields[2])
        self.nb_copies = float(fields[3])
        self.motif_size = float(fields[4])
        self.percNatches = float(fields[5])
        self.percIndels = float(fields[6])
        self.alignmentScore = float(fields[7])
        self.percA = float(fields[8])
        self.percC = float(fields[9])
        self.percG = float(fields[10])
        self.percT = float(fields[11])
        self.entropy = float(fields[12])
        self.motifSeq = str(fields[13])
        self.trSeq = str(fields[14])
        self.noTr5 = str(fields[15])
        self.noTr3 = str(fields[16])

        ## Not available by default
        self.trLen = self.end - self.beg
        self.seqLen = None

    def fraction_input_seq(self):
        '''
        Compute the fraction of the input sequence corresponding to the tandem repeat call
        '''
        # A) Compute fraction
        if self.seqLen is not None:
            trFraction = self.trLen / self.seqLen * 100

        # B) Not compute as length not available
        else:
            trFraction = None

        return trFraction
