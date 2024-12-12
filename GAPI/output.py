## DEPENDENCIES ##
# External
import pybedtools
import mappy as mp
import itertools
import statistics
from collections import Counter

# Internal
from GAPI import structures
from GAPI import formats
from GAPI import bamtools


def INS2VCF_junction(metaclusters, index, refLengths, source, build, species, outName, outDir):
    '''
    Write INS calls into a VCF file

    Input:
        1. metaclusters: list containing list of INS metaclusters
        2. index: minimap2 index for the reference genome 
        3. refLengths: Dictionary containing reference ids as keys and as values the length for each reference
        4. source: software version used to generate the insertion calls
        5. build: reference genome build
        6. species: specie
        7. outName: Output file name
        8. outDir: Output directory

    Output: vcf file containing identified metaclusters
    '''
    ## 1. Initialize VCF 
    VCF = formats.VCF()

    ## 2. Create header
    ## Define info 
    info = {'VTYPE': ['.', 'String', 'Type of variant'], \
            'ITYPE': ['.', 'String', 'Type of structural variant'], \
            'MECHANISM': ['.', 'String', 'Insertion mechanism'], \
            'FAM': ['.', 'String', 'Repeat family'], \
            'SUBFAM': ['.', 'String', 'Repeat subfamily'], \
            'GERMDB': ['.', 'String', 'List of germline variation databases where the variant is reported'], \
            'CIPOS': ['2', 'Integer', 'Confidence interval around POS for imprecise variants'], \
            'CYTOID': ['.', 'String', 'Source element cytoband identifier'], \
            'NBEXONS': ['1', 'Integer', 'Number of exons for a processed pseudogene insertion'], \
            'SRCGENE': ['.', 'String', 'Source gene for a processed psendogene insertion'], \
            'STRAND': ['.', 'String', 'Insertion DNA strand (+ or -)'], \
            'REGION': ['.', 'String', 'Genomic region where insertion occurs'], \
            'GENE': ['.', 'String', 'HUGO gene symbol'], \
            'REP': ['.', 'String', 'Families for annotated repeats at the insertion region'], \
            'REPSUB': ['.', 'String', 'Subfamilies for annotated repeats at the insertion region'], \
            'DIST': ['.', 'Integer', 'Distance between insertion breakpoint and annotated repeat'], \
            'NBTOTAL': ['1', 'Integer', 'Total number of insertion supporting reads'], \
            'NBTUMOR': ['1', 'Integer', 'Number of insertion supporting reads in the tumour'], \
            'NBNORMAL': ['1', 'Integer', 'Number of insertion supporting reads in the normal'], \
            'NBSPAN': ['1', 'Integer', 'Number of spanning supporting reads'], \
            'NBCLIP': ['1', 'Integer', 'Number of clipping supporting reads'], \
            'LEN': ['1', 'Integer', 'Insertion length'], \
            'CV': ['1', 'Float', 'Length coefficient of variation'], \
            'RTLEN': ['1', 'Integer', 'Inserted retrotransposon length'], \
            'TRUN5LEN': ['1', 'Integer', 'Size of 5prime truncation'], \
            'TRUN3LEN': ['1', 'Integer', 'Size of 3prime truncation'], \
            'FULL': ['0', 'Flag', 'Full length mobile element'], \
            'TDLEN': ['1', 'Integer', 'Transduction length'], \
            'INVLEN': ['1', 'Integer', '5-inversion length'], \
            'PERCR': ['1', 'Float', 'Percentage of inserted sequence that has been resolved'], \
            'QHITS': ['.', 'String', 'Coordinates for inserted sequence hits on the reference'], \
            'THITS': ['.', 'String', 'Inserted sequence hits on the reference'], \
            'RTCOORD': ['.', 'String', 'Coordinates for inserted retrotransposon piece of sequence'], \
            'POLYA': ['0', 'Flag', 'PolyA tail identified'], \
            'INSEQ': ['.', 'String', 'Inserted sequence'], \
            }
            
    ## Create header
    VCF.create_header(source, build, species, refLengths, info, {}, [])

    ## 3. Add insertion calls to the VCF
    ## 3.1 Load reference index
    reference = mp.Aligner(fn_idx_in=index) # comment as time consuming

    ## 3.2 Iterate over INS metaclusters
    for metacluster in metaclusters:

        print ('metacluster.bridge ' + str(metacluster.bridge))

        ## Collect insertion basic features
        CHROM = metacluster.ref
        POS, CIPOS = metacluster.mean_pos()
        ID = '.'
        REF = reference.seq(CHROM, POS, POS + 1)
        ALT = '<INS>'
        QUAL = '.'
        FILTER = 'PASS' if not metacluster.failedFilters else ','.join(metacluster.failedFilters)
        
        ## Collect extra insertion features to include at info field
        INFO = {}
        repeats = metacluster.repeatAnnot if hasattr(metacluster, 'repeatAnnot') else []        

        INFO['VTYPE'] = metacluster.mutOrigin
        INFO['ITYPE'] = metacluster.SV_features['INS_TYPE'] if 'INS_TYPE' in metacluster.SV_features else None
        INFO['MECHANISM'] = metacluster.SV_features['MECHANISM'] if 'MECHANISM' in metacluster.SV_features else None        
        INFO['FAM'] = ','.join(metacluster.SV_features['FAMILY']) if ('FAMILY' in metacluster.SV_features and metacluster.SV_features['FAMILY']) else None
        INFO['SUBFAM'] = ','.join(metacluster.SV_features['SUBFAMILY']) if ('SUBFAMILY' in metacluster.SV_features and metacluster.SV_features['SUBFAMILY']) else None
        INFO['GERMDB'] = metacluster.germlineDb if hasattr(metacluster, 'germlineDb') else None       
        INFO['CIPOS'] = str(CIPOS[0]) + ',' + str(CIPOS[1]) 
        INFO['CYTOID'] = ','.join(metacluster.SV_features['CYTOBAND']) if ('CYTOBAND' in metacluster.SV_features and metacluster.SV_features['CYTOBAND']) else None
        INFO['NBEXONS'] = metacluster.SV_features['NB_EXONS'] if 'NB_EXONS' in metacluster.SV_features else None
        INFO['SRCGENE'] = ','.join(metacluster.SV_features['SOURCE_GENE']) if 'SOURCE_GENE' in metacluster.SV_features else None
        INFO['STRAND'] = metacluster.SV_features['STRAND'] if 'STRAND' in metacluster.SV_features else None
        INFO['REGION'], INFO['GENE'] = metacluster.geneAnnot if hasattr(metacluster, 'geneAnnot') else (None, None)
        INFO['REP'] = ','.join([repeat['family'] for repeat in repeats]) if repeats else None 
        INFO['REPSUB'] = ','.join([repeat['subfamily'] for repeat in repeats]) if repeats else None   
        INFO['DIST'] = ','.join([str(repeat['distance']) for repeat in repeats]) if repeats else None
        INFO['NBTOTAL'], INFO['NBTUMOR'], INFO['NBNORMAL'] = str(metacluster.nbTotal), str(metacluster.nbTumour), str(metacluster.nbNormal) 
        INFO['NBSPAN'], INFO['NBCLIP'] = str(metacluster.nbINS), str(metacluster.nbCLIPPING)
        INFO['LEN'] = metacluster.consensusEvent.length if metacluster.consensusEvent is not None else None
        INFO['CV'] = metacluster.cv
        INFO['RTLEN'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        INFO['TRUN5LEN'] = metacluster.SV_features['TRUNCATION_5_LEN'] if 'TRUNCATION_5_LEN' in metacluster.SV_features else None
        INFO['TRUN3LEN'] = metacluster.SV_features['TRUNCATION_3_LEN'] if 'TRUNCATION_3_LEN' in metacluster.SV_features else None
        INFO['FULL'] = metacluster.SV_features['IS_FULL'] if 'IS_FULL' in metacluster.SV_features else None
        INFO['TDLEN'] = metacluster.SV_features['TRANSDUCTION_LEN'] if 'TRANSDUCTION_LEN' in metacluster.SV_features else None
        INFO['INVLEN'] = metacluster.SV_features['INVERSION_LEN'] if 'INVERSION_LEN' in metacluster.SV_features else None
        INFO['PERCR'] = metacluster.SV_features['PERC_RESOLVED'] if 'PERC_RESOLVED' in metacluster.SV_features else None
        INFO['QHITS'] = None if metacluster.insertHits is None else ','.join([ 'insertedSeq' + ':' + str(alignment.qBeg) + '-' + str(alignment.qEnd) for alignment in metacluster.insertHits.alignments ])
        INFO['THITS'] = None if metacluster.insertHits is None else ','.join([ alignment.tName + ':' + str(alignment.tBeg) + '-' + str(alignment.tEnd) for alignment in metacluster.insertHits.alignments ])
        INFO['RTCOORD'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        INFO['POLYA'] = metacluster.SV_features['POLYA'] if 'POLYA' in metacluster.SV_features else None
        #INFO['INSEQ'] = metacluster.consensusEvent.pick_insert() if metacluster.consensusEvent is not None else None

        ## Create VCF variant object
        fields = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, {}]

        ## Add variant to the VCF
        INS = formats.VCF_variant(fields)
        VCF.add(INS)
        
    ## 4. Sort VCF
    VCF.sort()

    ## 5. Write VCF in disk
    IDS = ['VTYPE', 'ITYPE', 'MECHANISM', 'FAM', 'SUBFAM', 'GERMDB', 'CIPOS', 'CYTOID', \
           'NBEXONS', 'SRCGENE', 'STRAND', 'REGION', 'GENE', 'REP', 'REPSUB', 'DIST', \
           'NBTOTAL', 'NBTUMOR', 'NBNORMAL', 'NBSPAN', 'NBCLIP', 'LEN', 'CV', 'RTLEN', \
           'TRUN5LEN', 'TRUN3LEN', 'FULL', 'TDLEN', 'INVLEN', 'PERCR', \
           'QHITS', 'THITS', 'RTCOORD', 'POLYA', 'INSEQ']

    VCF.write(IDS, [], outName, outDir)


def INS2VCF(metaclusters, index, refLengths, source, build, species, outName, outDir):
    '''
    Write INS calls into a VCF file

    Input:
        1. metaclusters: list containing list of INS metaclusters
        2. index: minimap2 index for the reference genome 
        3. refLengths: Dictionary containing reference ids as keys and as values the length for each reference
        4. source: software version used to generate the insertion calls
        5. build: reference genome build
        6. species: specie
        7. outName: Output file name
        8. outDir: Output directory

    Output: vcf file containing identified metaclusters
    '''
    ## 1. Initialize VCF 
    VCF = formats.VCF()

    ## 2. Create header
    ## Define info 
    info = {'VTYPE': ['.', 'String', 'Type of variant'], \
            'ITYPE': ['.', 'String', 'Type of structural variant'], \
            'MECHANISM': ['.', 'String', 'Insertion mechanism'], \
            'FAM': ['.', 'String', 'Repeat family'], \
            'SUBFAM': ['.', 'String', 'Repeat subfamily'], \
            'GERMDB': ['.', 'String', 'List of germline variation databases where the variant is reported'], \
            'CIPOS': ['2', 'Integer', 'Confidence interval around POS for imprecise variants'], \
            'CYTOID': ['.', 'String', 'Source element cytoband identifier'], \
            'NBEXONS': ['1', 'Integer', 'Number of exons for a processed pseudogene insertion'], \
            'SRCGENE': ['.', 'String', 'Source gene for a processed psendogene insertion'], \
            'STRAND': ['.', 'String', 'Insertion DNA strand (+ or -)'], \
            'REGION': ['.', 'String', 'Genomic region where insertion occurs'], \
            'GENE': ['.', 'String', 'HUGO gene symbol'], \
            'REP': ['.', 'String', 'Families for annotated repeats at the insertion region'], \
            'REPSUB': ['.', 'String', 'Subfamilies for annotated repeats at the insertion region'], \
            'DIST': ['.', 'Integer', 'Distance between insertion breakpoint and annotated repeat'], \
            'NBTOTAL': ['1', 'Integer', 'Total number of insertion supporting reads'], \
            'NBTUMOR': ['1', 'Integer', 'Number of insertion supporting reads in the tumour'], \
            'NBNORMAL': ['1', 'Integer', 'Number of insertion supporting reads in the normal'], \
            'NBSPAN': ['1', 'Integer', 'Number of spanning supporting reads'], \
            'NBCLIP': ['1', 'Integer', 'Number of clipping supporting reads'], \
            'LEN': ['1', 'Integer', 'Insertion length'], \
            'CV': ['1', 'Float', 'Length coefficient of variation'], \
            'RTLEN': ['1', 'Integer', 'Inserted retrotransposon length'], \
            'TRUN5LEN': ['1', 'Integer', 'Size of 5prime truncation'], \
            'TRUN3LEN': ['1', 'Integer', 'Size of 3prime truncation'], \
            'FULL': ['0', 'Flag', 'Full length mobile element'], \
            'TDLEN': ['1', 'Integer', 'Transduction length'], \
            'INVLEN': ['1', 'Integer', '5-inversion length'], \
            'PERCR': ['1', 'Float', 'Percentage of inserted sequence that has been resolved'], \
            'QHITS': ['.', 'String', 'Coordinates for inserted sequence hits on the reference'], \
            'THITS': ['.', 'String', 'Inserted sequence hits on the reference'], \
            'RTCOORD': ['.', 'String', 'Coordinates for inserted retrotransposon piece of sequence'], \
            'POLYA': ['0', 'Flag', 'PolyA tail identified'], \
            'INSEQ': ['.', 'String', 'Inserted sequence'], \
            }
            
    ## Create header
    VCF.create_header(source, build, species, refLengths, info, {}, None)

    ## 3. Add insertion calls to the VCF
    ## 3.1 Load reference index
    reference = mp.Aligner(fn_idx_in=index) # comment as time consuming

    ## 3.2 Iterate over INS metaclusters
    for metacluster in metaclusters:

        ## Collect insertion basic features
        CHROM = metacluster.ref
        POS, CIPOS = metacluster.mean_pos()
        ID = '.'
        REF = reference.seq(CHROM, POS, POS + 1)
        ALT = '<INS>'
        QUAL = '.'
        FILTER = 'PASS' if not metacluster.failedFilters else ','.join(metacluster.failedFilters)
        
        ## Collect extra insertion features to include at info field
        INFO = {}
        repeats = metacluster.repeatAnnot if hasattr(metacluster, 'repeatAnnot') else []        

        INFO['VTYPE'] = metacluster.mutOrigin
        INFO['ITYPE'] = metacluster.SV_features['INS_TYPE'] if 'INS_TYPE' in metacluster.SV_features else None
        INFO['MECHANISM'] = metacluster.SV_features['MECHANISM'] if 'MECHANISM' in metacluster.SV_features else None        
        INFO['FAM'] = ','.join(metacluster.SV_features['FAMILY']) if ('FAMILY' in metacluster.SV_features and metacluster.SV_features['FAMILY']) else None
        INFO['SUBFAM'] = ','.join(metacluster.SV_features['SUBFAMILY']) if ('SUBFAMILY' in metacluster.SV_features and metacluster.SV_features['SUBFAMILY']) else None
        INFO['GERMDB'] = metacluster.germlineDb if hasattr(metacluster, 'germlineDb') else None       
        INFO['CIPOS'] = str(CIPOS[0]) + ',' + str(CIPOS[1]) 
        INFO['CYTOID'] = ','.join(metacluster.SV_features['CYTOBAND']) if ('CYTOBAND' in metacluster.SV_features and metacluster.SV_features['CYTOBAND']) else None
        INFO['NBEXONS'] = metacluster.SV_features['NB_EXONS'] if 'NB_EXONS' in metacluster.SV_features else None
        INFO['SRCGENE'] = ','.join(metacluster.SV_features['SOURCE_GENE']) if 'SOURCE_GENE' in metacluster.SV_features else None
        INFO['STRAND'] = metacluster.SV_features['STRAND'] if 'STRAND' in metacluster.SV_features else None
        INFO['REGION'], INFO['GENE'] = metacluster.geneAnnot if hasattr(metacluster, 'geneAnnot') else (None, None)
        INFO['REP'] = ','.join([repeat['family'] for repeat in repeats]) if repeats else None 
        INFO['REPSUB'] = ','.join([repeat['subfamily'] for repeat in repeats]) if repeats else None   
        INFO['DIST'] = ','.join([str(repeat['distance']) for repeat in repeats]) if repeats else None
        INFO['NBTOTAL'], INFO['NBTUMOR'], INFO['NBNORMAL'] = str(metacluster.nbTotal), str(metacluster.nbTumour), str(metacluster.nbNormal) 
        INFO['NBSPAN'], INFO['NBCLIP'] = str(metacluster.nbINS), str(metacluster.nbCLIPPING)
        INFO['LEN'] = metacluster.consensusEvent.length if metacluster.consensusEvent is not None else None
        INFO['CV'] = metacluster.cv
        INFO['RTLEN'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        INFO['TRUN5LEN'] = metacluster.SV_features['TRUNCATION_5_LEN'] if 'TRUNCATION_5_LEN' in metacluster.SV_features else None
        INFO['TRUN3LEN'] = metacluster.SV_features['TRUNCATION_3_LEN'] if 'TRUNCATION_3_LEN' in metacluster.SV_features else None
        INFO['FULL'] = metacluster.SV_features['IS_FULL'] if 'IS_FULL' in metacluster.SV_features else None
        INFO['TDLEN'] = metacluster.SV_features['TRANSDUCTION_LEN'] if 'TRANSDUCTION_LEN' in metacluster.SV_features else None
        INFO['INVLEN'] = metacluster.SV_features['INVERSION_LEN'] if 'INVERSION_LEN' in metacluster.SV_features else None
        INFO['PERCR'] = metacluster.SV_features['PERC_RESOLVED'] if 'PERC_RESOLVED' in metacluster.SV_features else None
        INFO['QHITS'] = None if metacluster.insertHits is None else ','.join([ 'insertedSeq' + ':' + str(alignment.qBeg) + '-' + str(alignment.qEnd) for alignment in metacluster.insertHits.alignments ])
        INFO['THITS'] = None if metacluster.insertHits is None else ','.join([ alignment.tName + ':' + str(alignment.tBeg) + '-' + str(alignment.tEnd) for alignment in metacluster.insertHits.alignments ])
        INFO['RTCOORD'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        INFO['POLYA'] = metacluster.SV_features['POLYA'] if 'POLYA' in metacluster.SV_features else None
        INFO['INSEQ'] = metacluster.consensusEvent.pick_insert() if metacluster.consensusEvent is not None else None

        ## Create VCF variant object
        fields = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, {}]

        ## Add variant to the VCF
        INS = formats.VCF_variant(fields)
        VCF.add(INS)
        
    ## 4. Sort VCF
    VCF.sort()

    ## 5. Write VCF in disk
    infoIds = ['VTYPE', 'ITYPE', 'MECHANISM', 'FAM', 'SUBFAM', 'GERMDB', 'CIPOS', 'CYTOID', \
           'NBEXONS', 'SRCGENE', 'STRAND', 'REGION', 'GENE', 'REP', 'REPSUB', 'DIST', \
           'NBTOTAL', 'NBTUMOR', 'NBNORMAL', 'NBSPAN', 'NBCLIP', 'LEN', 'CV', 'RTLEN', \
           'TRUN5LEN', 'TRUN3LEN', 'FULL', 'TDLEN', 'INVLEN', 'PERCR', \
           'QHITS', 'THITS', 'RTCOORD', 'POLYA', 'INSEQ']

    VCF.write(infoIds, [], outName, outDir)

def alt_assemblies2fasta(metaclustersPass, metaclustersFailed, offset, outDir):
    '''
    Write fasta file with the consensus assemblies supporting the variant
    
    Input:
        1. metaclustersPass: List of metaclusters that passed all filters
        2. metaclustersFailed: List of metaclusters that do not passed all filters
        3. offset: Number of nt surrounding the event
        4. outDir: Output directory
    
    Output:
        Fasta file (alt-assembly.fa) in outDir 
    '''
    # Set variables and create fasta object. 
    offset = 1000
    altAssemblies_fa = outDir + '/' + 'alt-assembly.fa'
    altAssemblies = formats.FASTA()
    
    # Store all info in fasta object seqDict 
    if 'INS' in metaclustersPass:
        
        for metacluster in metaclustersPass['INS']:
            
            # create header using ref_bkp_PASS as dict key
            svId = '_'.join([str(metacluster.ref), str(metacluster.mean_pos()[0]), 'PASS'])
            # keep sequence as dict value
            altAssemblies.seqDict[svId] = metacluster.consensusEvent.pick_ALT(offset) if metacluster.consensusEvent is not None else None

    if 'INS' in metaclustersFailed:
        
        for metacluster in metaclustersPass['INS']:
            
            # create header using ref_bkp_PASS as dict key
            svId = '_'.join([str(metacluster.ref), str(metacluster.mean_pos()[0]), 'FAILED'])
            # keep sequence as dict value
            altAssemblies.seqDict[svId] = metacluster.consensusEvent.pick_ALT(offset) if metacluster.consensusEvent is not None else None
    
    # Write fasta object dict to fasta file                
    altAssemblies.write(altAssemblies_fa)