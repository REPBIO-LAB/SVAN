'''
Module 'databases' - Contains functions and classes to deal with different types of databases
'''

## DEPENDENCIES ##
# External
import subprocess

# Internal
from GAPI import unix 
from GAPI import log

## FUNCTIONS ##
def buildRetrotransposonDb(fastaDir, includeTransduced, outDir):
    '''
    Build database containing retrotransposon related sequences (consensus sequences, transduced regions, ...)

    Input:
        1. fastaDir: Directory containing reference fasta files (retrotransposon consensus sequences, transduced regions...)
        2. includeTransduced: Boolean to specify if include (True) or not (False) transduced regions from known source elements in the database
        3. outDir: Output directory

    Output:
        1. retrotransposonDb: Fasta file containing retrotransposon related sequences 
    '''

    ## 0. Create logs directory
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Create database fasta file ##
    retrotransposonDb = outDir + '/retrotransposonDb.fa'

    # A) Include retrotransposon consensus sequences + transduced regions
    if includeTransduced:

        consensusDb = fastaDir + 'consensusDb.fa'
        transducedDb = fastaDir + 'transducedDb.masked.fa'
        files = [consensusDb, transducedDb]

        with open(retrotransposonDb, 'w') as outFile:
            # Iterate over each fasta and write fasta into output database
            for f in files:
                with open(f) as inFile:
                    outFile.write(inFile.read())
        
    # B) Only include retrotransposon consensus sequences
    else:
        consensusDb = fastaDir + 'consensusDb.fa'

        with open(retrotransposonDb, 'w') as outFile:
            with open(consensusDb) as inFile:
                outFile.write(inFile.read())

    ## 2. Index retrotransposon database fasta file ##
    index = outDir + '/retrotransposonDb.mmi'
    err = open(logDir + '/index.err', 'w') 
    command = 'minimap2 -k 10 -w 1 -d ' + index + ' ' + retrotransposonDb 
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'BUILD-DATABASE'
        msg = 'Database indexing failed' 
        log.step(step, msg)

    return retrotransposonDb, index


def buildIdentityDb(metacluster, db, outDir):

    # Coger la identity del primer evento DISCORDANT que aparezca (pq los clipping no tienen identity)
    identity = next(event.identity for event in metacluster.events if event.type == "DISCORDANT")
    
    specificIdentity = max(set([event.specificIdentity for event in metacluster.events if event.type == "DISCORDANT"]), key=[event.specificIdentity for event in metacluster.events if event.type == "DISCORDANT"].count)

    specificHeader = '"consensus|' + specificIdentity + '|' + identity + '"'

    #dbSpecificIdentity = outDir + '/' + str(CLIPPING_clusterID) + '_specificIdentity.fa'
    dbSpecificIdentity = outDir + '/' + str(metacluster.id) + '_specificIdentity.fa' 

    # Build database con identity + ref
    # TODO coger la identity que ya hemos reconocido, en lugar de toda la db!
    # 2. Cojo de la db de identities la secuencia que fue asignada como identity:
    # TODO poner bien el status y todo eso
    #err = open(logDir + '/index.err', 'w') 
    command = 'samtools faidx  ' + db + ' ' + specificHeader + ' -o ' + dbSpecificIdentity
    status = subprocess.call(command, shell=True)
    # indexo
    indexDbSpecificIdentity = dbSpecificIdentity.replace("_specificIdentity.fa", "_refIdentityDb.mmi")

    ## DESILENCIAAAAR
    # TODO
    # ponerlo bien
    #err = open(logDir + '/index.err', 'w')
    
    command = 'minimap2 -k 10 -w 1 -d ' + indexDbSpecificIdentity + ' ' + dbSpecificIdentity 
    status = subprocess.call(command, shell=True)

    if status != 0:
        step = 'BUILD-VIRUS-DATABASE'
        msg = 'Database indexing failed' 
        log.step(step, msg)
        
    return indexDbSpecificIdentity

def create_transduced_bed(sourceBed, srcEnds, size, buffer, outDir):
    '''
    Create bed file containing regions frequently transduced by source elements
    
    Input:
        1. sourceBed: Bed file source elements coordinates, cytoband identifier, family and orientation. Following fields required:
                      1) chrom
                      2) beg
                      3) end
                      4) cytobandId
                      5) family
                      6) strand
        2. srcEnds: source elements ends to look for transductions. [3, 5]
        3. size: transduced region size
        4. buffer: buffer to apply to the end of the elemnt. ME end - buffer to define transduced region beg
        5. outDir: Output directory
        
    Output:
        1. transducedPath: Bed file containing transduced region coordinates
    '''
    
    ## Open file handlers
    sourceBed = open(sourceBed, 'r')
    transducedPath = outDir + '/transduced_regions.bed'
    transducedBed = open(transducedPath, 'w')
    
    ## Write header in outfile
    row = "\t".join(['#ref', 'tdBeg', 'tdEnd', 'name', 'cytobandId', 'family', 'strand', 'srcEnd', "\n"])
    transducedBed.write(row)
    
    ## Read bed with source elements annotation line by line
    for line in sourceBed:
        line = line.rstrip('\r\n')
        
        ## Discard header
        if not line.startswith("#"):   
            
            fieldsList = line.split("\t")
            ref, beg, end, name, cytobandId, family, strand = fieldsList
                    
            ## For each end of the source elements
            for srcEnd in srcEnds:
                
                ## a) Elements:
                # beg ---------------> end ........transduced........ end + size
                # beg <--------------- end ........transduced........ end + size
                if ((strand == '+') and (srcEnd == '3')) or ((strand == '-') and (srcEnd == '5')):
                    
                    tdBeg = int(end) - buffer
                    tdEnd = int(end) + size
                    
                ## b) Elements:
                # beg - size ........transduced........ beg <--------------- end
                # beg - size ........transduced........ beg ---------------> end
                elif ((strand == '-') and (srcEnd == '3')) or ((strand == '+') and (srcEnd == '5')):
                    
                    tdBeg = int(beg) - size
                    tdEnd = int(beg) + buffer
                    
                else:
                    print("There is a problem with the transduction coordinates")
                    
                ## Write into output file
                row = "\t".join([ref, str(tdBeg), str(tdEnd), name, cytobandId, family, strand, str(srcEnd), "\n"])
                transducedBed.write(row)
            
    return transducedPath