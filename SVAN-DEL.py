## DEPENDENCIES ##
# External
import os
import sys
import argparse
import copy
import re
import pandas as pd
import pysam
from collections import Counter
import os.path
import subprocess
import os.path

# Internal
import formats
import alignment
import unix
import sequences
import gRanges
import retrotransposons
import structures
import log

###############
## Functions ##
###############
def annot_vntr(VCF, TRF_out, vntrBinDb):
    '''
    '''
    ### Load TRF calls
    TRF = formats.TRF()
    TRF.read(TRF_out)

    ## Create fasta object containing insertion sequences
    FASTA_ins = ins2fasta(VCF)

    ## Add sequence size to the TRF objects
    TRF.add_seq_size(FASTA_ins)

    ## Create VCF with VNTR insertion calls
    vntr_VCF = formats.VCF()
    vntr_VCF.header = VCF.header

    vntr_VCF.info_order = VCF.info_order
    vntr_VCF.format_order = VCF.format_order
    
    ## Create VCF with non-VNTR insertion
    NO_vntr_VCF = formats.VCF()
    NO_vntr_VCF.header = VCF.header

    NO_vntr_VCF.info_order = VCF.info_order
    NO_vntr_VCF.format_order = VCF.format_order

    ## Assess for each variant in the VCF file if VNTR
    for variant in VCF.variants:
 
        insId = variant.ID + '_' + variant.chrom + ':' + str(variant.pos) 

        ## Check if insertion overlap a VNTR locus
        #buffer = 10 
        #overlaps = vntrBinDb[variant.chrom].collect_interval(variant.pos - buffer, variant.pos + buffer, ['VNTR'])

        #if overlaps: 
        #    tmp = {'IN_VNTR_DB': True}
        #else:
        #    tmp = {'IN_VNTR_DB': False}

        #variant.info.update(tmp)

        ## Analysis of sequence composition
        if insId not in TRF.callsDict:
            vntr_bool = False
            perc = 0
        else:
            vntr_bool, perc, nbMotifs, motifs = TRF.callsDict[insId].is_tandem_repeat(75, 10, 80)

        ## A) Call VNTR
        if vntr_bool:

            ## Add VNTR annotation to the variant info
            tmp = {'DTYPE_N': 'VNTR', 'PERC_RESOLVED': perc, 'NB_MOTIFS': nbMotifs, 'MOTIFS': motifs}
            variant.info.update(tmp)

            ## Add variant to the VCF with VNTR
            vntr_VCF.add(variant)

        ## B) NO VNTR
        else:
            NO_vntr_VCF.add(variant)

    return vntr_VCF, NO_vntr_VCF    

def align2flankingRegion(variant, refIndex, offset, outDir):
    '''
    '''
    ## Insertion identifier
    insId = variant.ID + '_' + variant.chrom + ':' + str(variant.pos) 
    
    ## Create fasta with the insert
    fasta_insert = formats.FASTA() 
    seq_insert = variant.ref
    fasta_insert.seqDict['insert' + '_' + insId] = seq_insert
    insert_path = outDir + '/insert.fa'
    fasta_insert.write(insert_path)

    ## Collect sequence flanking the variant breakpoint
    insLen = len(variant.ref)
    beg = variant.pos - insLen - offset
    end = variant.pos + insLen + offset

    if beg < 1:
        beg = 1

    coord = variant.chrom + ':' + str(beg) + '-' + str(end)
    seq_flank = refIndex.fetch(region=coord)

    ## Create fasta with breakpoint flanking sequence
    fasta_flank = formats.FASTA() 
    fasta_flank.seqDict['flank' + '_' + insId] = seq_flank
    flank_path = outDir + '/flanking_region.fa'

    if not os.path.isfile(flank_path):

        fasta_flank.write(flank_path)
        alignment.index_bwa(flank_path, outDir)

    ## Align
    PAF_path = outDir + '/insert2flankingRegion.paf'

    if not os.path.isfile(PAF_path):
        SAM_path = alignment.alignment_bwa(insert_path, flank_path, 'insert2flankingRegion', 2, outDir)     
        PAF_path = alignment.sam2paf(SAM_path, 'insert2flankingRegion', outDir)

    ## Read PAF
    PAF_region = formats.PAF()

    if os.path.isfile(PAF_path):
        PAF_region.read(PAF_path)

    ## Cleanup
    #unix.rm([insert_path, flank_path, flank_path + '.bwt', flank_path + '.pac', flank_path + '.ann', flank_path + '.amb', flank_path + '.sa', SAM_path, outDir + '/index.err', outDir + '/align.err', outDir + '/sam2paf.err'])

    return PAF_region

def annot_retro(VCF, exonsBinDb, PAFs_ref, reference, repeatsDb, TRF_out, outDir):
    '''
    '''
    ## 1. Identify RETRO candidates 
    ################################
    # == INS with a poly(A/T) at one of its ends
    candidate_VCF, NO_RETRO_VCF, retro_annot, stats = call_retro_candidate(VCF)

    ## 2. Report some stats 
    ########################
    ## Classification
    for i in ['RETRO_CANDIDATE', 'HOMO_A', 'HOMO_T', 'NO_TAIL', 'DOUBLE_TAIL']:
        print('NB_' + i, ': ', stats[i])

    ## Strand for RETRO candidates
    counts = pd.Series([retro_annot[i]['STRAND'] for i in retro_annot.keys()]).value_counts()
    print('STRAND_CANDIDATES: ', counts)

    ## 3. Annotate and trim SVA-like hexamer tandem repeats at the 5' end of insert
    #################################################################################
    annotate_hexamer(candidate_VCF, retro_annot, TRF_out)

    ## 4. Annotate and trim TSD for candidate RETRO
    ##############################################
    annotate_tsd(candidate_VCF, reference, retro_annot)

    ## 5. Annotate Endonuclease cleavage site for candidate RETRO
    ##############################################################
    #annotate_endomotif(candidate_VCF, reference, retro_annot)

    ## 6. Reverse candidate RETRO sequences in - orientation
    #######################################################
    ## Note: work with sequence after TSD trimming 
    rev2forward(retro_annot)

    ## 7. Processed pseudogene annotation
    ######################################
    for insId in retro_annot:

        PAF_ref = PAFs_ref[insId] if insId in PAFs_ref else formats.PAF()

        retro_annot[insId] = processed_pseudogene_annot(retro_annot[insId], PAF_ref, exonsBinDb)
    
    ## 8. Search for 3-prime transductions and trim polyA + transduced sequences
    #############################################################################
    search4transductions_3prime(retro_annot, exonsBinDb, reference, outDir)

    ## 9. Annot retrotransposed element
    #####################################
    retro_annot_pass, retro_annot_fail = annot_retroelement(retro_annot, exonsBinDb, repeatsDb, TRF_out, reference, outDir)
        
    print('NB_RETRO_PASS: ', len(retro_annot_pass))
    print('NB_RETRO_FAIL: ', len(retro_annot_fail))

    '''
    for insId in retro_annot_pass:
        if '5PRIME_MOB_GENE' in retro_annot_pass[insId]:
            print('PASS_5PRIME_MOB_GENE: ', retro_annot_pass[insId]['5PRIME_MOB_GENE'])
        
        if '3PRIME_MOB_GENE' in retro_annot_pass[insId]:
            print('PASS_3PRIME_MOB_GENE: ', retro_annot_pass[insId]['3PRIME_MOB_GENE'])

    for insId in retro_annot_fail:
        if '5PRIME_MOB_GENE' in retro_annot_fail[insId]:
            print('FAIL_5PRIME_MOB_GENE: ', retro_annot_fail[insId]['5PRIME_MOB_GENE'])
        
        if '3PRIME_MOB_GENE' in retro_annot_fail[insId]:
            print('FAIL_3PRIME_MOB_GENE: ', retro_annot_fail[insId]['3PRIME_MOB_GENE'])
    '''
          
    ## 10. Separate insertions into RETRO and NO_RETRO vcf files
    ########################################################
    ## Create VCF for RETRO insertion calls
    RETRO_VCF = formats.VCF()
    RETRO_VCF.header = NO_RETRO_VCF.header
    RETRO_VCF.info_order = NO_RETRO_VCF.info_order
    RETRO_VCF.format_order = NO_RETRO_VCF.format_order

    ## Iterate over variants in the candidate VCF assigning them to the RETRO and no_RETRO VCFs
    for variant in candidate_VCF.variants:
 
        insId = variant.ID + '_' + variant.chrom + ':' + str(variant.pos) 

        ## A) NO RETRO
        if insId in retro_annot_fail:

            NO_RETRO_VCF.add(variant)

        ## B) RETRO
        if insId in retro_annot_pass:

            ## Add annotation to the VCF info field
            variant.info.update(retro_annot_pass[insId])

            ## Add RETRO to the VCF
            RETRO_VCF.add(variant)

    return RETRO_VCF, NO_RETRO_VCF

def annot_mei_not_retro(VCF, PAFs_repeatDb):
    '''
    Include structure annotation
    '''
    ## Create VCF with not retrotransposed MEI
    MEI_NOT_RETRO_VCF = formats.VCF()
    MEI_NOT_RETRO_VCF.header = VCF.header

    MEI_NOT_RETRO_VCF.info_order = VCF.info_order
    MEI_NOT_RETRO_VCF.format_order = VCF.format_order
    
    ## Create VCF without MEI
    NO_MEI_VCF = formats.VCF()
    NO_MEI_VCF.header = VCF.header

    NO_MEI_VCF.info_order = VCF.info_order
    NO_MEI_VCF.format_order = VCF.format_order

    ## For each variant assess if retrotransposon insert without retrotransposition hallmarks
    for variant in VCF.variants:
 
        insId = variant.ID + '_' + variant.chrom + ':' + str(variant.pos) 

        ## Discard if no hits on retrotransposon templates 
        if (insId not in PAFs_repeatDb) :
            NO_MEI_VCF.add(variant)
            continue

        ## Copy PAF
        hits_PAF = copy.deepcopy(PAFs_repeatDb[insId])

        ## Filter hits based on mapping quality
        hits_PAF.alignments = [hit for hit in hits_PAF.alignments if hit.MAPQ >= 20]

        ## Discard if no hit passes the filters
        if not hits_PAF.alignments:
            NO_MEI_VCF.add(variant)
            continue

        ## Create alignment chain
        chain = hits_PAF.chain(10, 25)

        ## Call not retrotransposed MEI if:
        # - a) Chain covers at least 75% of the query sequence
        if (chain.perc_query_covered() >= 75):

            ## Collect info regarding the templates in the chain
            nbTemplates, templates = chain.nb_templates()
            templates = [i.split('|')[1] for i in templates]

            ## Annotate structure and conformation for Alu and L1
            if ('L1' in templates) or ('Alu' in templates):
                
                annot = {}
                annot['CHAIN_repeats'] = chain
                beg, end =  annot['CHAIN_repeats'].interval()
                annot['RT_LEN'] = end - beg
                annot = annotate_structure(annot)
                variant.info.update(annot)
            
            templates = ','.join(templates)
            tmp = {'DTYPE_N': 'solo_no_canonical', 'FAM_N': templates, 'PERC_RESOLVED': chain.perc_query_covered()}
            variant.info.update(tmp)
            
            ### Add call to the VCF
            MEI_NOT_RETRO_VCF.add(variant)

        else:
            NO_MEI_VCF.add(variant)

    return MEI_NOT_RETRO_VCF, NO_MEI_VCF


def call_retro_candidate(VCF):
    '''
    '''
    ## Create VCF with candidate RETRO calls
    candidate_VCF = formats.VCF()
    candidate_VCF.header = VCF.header
    candidate_VCF.info_order = VCF.info_order
    candidate_VCF.format_order = VCF.format_order
  
    ## Create VCF with non-candidate RETRO
    no_candidate_VCF = formats.VCF()
    no_candidate_VCF.header = VCF.header
    no_candidate_VCF.info_order = VCF.info_order
    no_candidate_VCF.format_order = VCF.format_order

    ## Create RETRO annot dictionary
    retro_annot = {}

    ## Create stats dict
    stats = {}
    stats['NO_TAIL'] = 0
    stats['DOUBLE_TAIL'] = 0
    stats['HOMO_A'] = 0
    stats['HOMO_T'] = 0
    stats['RETRO_CANDIDATE'] = 0

    ## For each variant
    for variant in VCF.variants:

        insId = variant.ID + '_' + variant.chrom + ':' + str(variant.pos) 

        ## Search for poly(A) tail at inserted sequence
        polyA, monomerA = search4polyA(variant.ref)

        ## Search for poly(T) tail at inserted sequence
        polyT, monomerT = search4polyT(variant.ref)

        ## A) Filter calls if poly(A) nor poly(T) found
        if not polyA and not polyT:
            stats['NO_TAIL'] += 1
            no_candidate_VCF.add(variant)
        
        ## B) Both poly(A) and poly(T) found
        elif polyA and polyT:
            stats['DOUBLE_TAIL'] += 1      
            no_candidate_VCF.add(variant)

        ## C) Poly(A) found
        elif polyA:
            percLen =  monomerA.length() / len(variant.ref) * 100
            noPolyLen = len(variant.ref) - monomerA.length()

            # a) A homopolymer
            if (percLen >= 95) or (noPolyLen < 10):
                stats['HOMO_A'] += 1            
                no_candidate_VCF.add(variant)

            # b) RETRO candidate
            else:
                retro_annot[insId] = {}
                retro_annot[insId]['DTYPE_N'] = 'UNK'
                retro_annot[insId]['STRAND'] = '+'
                retro_annot[insId]['DEL_LEN'] = len(variant.ref)
                retro_annot[insId]['POLY_MONOMER'] = monomerA
                stats['RETRO_CANDIDATE'] += 1   
                candidate_VCF.add(variant)
                
        ## D) Poly(T) found
        elif polyT:
            percLen =  monomerT.length() / len(variant.ref) * 100
            noPolyLen = len(variant.ref) - monomerT.length()

            # a) T homopolymer
            if (percLen >= 95) or (noPolyLen < 10):
                stats['HOMO_T'] += 1
                no_candidate_VCF.add(variant)

            # b) RETRO candidate
            else:
                retro_annot[insId] = {}
                retro_annot[insId]['DTYPE_N'] = 'UNK'
                retro_annot[insId]['STRAND'] = '-'  
                retro_annot[insId]['DEL_LEN'] = len(variant.ref)
                retro_annot[insId]['POLY_MONOMER'] = monomerT
                stats['RETRO_CANDIDATE'] += 1
                candidate_VCF.add(variant)

    return candidate_VCF, no_candidate_VCF, retro_annot, stats

def search4polyA(sequence):
    '''
    Search for poly(A) at sequence end 
    '''
    ## Configuration for monomere search:
    windowSize = 10
    maxWindowDist = 1
    minMonomerSize = 10
    minPurity = 80

    ## Seach poly(A) monomers
    targetMonomer = 'A'
    monomersA = sequences.find_monomers_backward(sequence, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

    if not monomersA:
        return False, None
        
    ## Select monomer closest to sequence end
    candidate = monomersA[-1]

    ## Filter out monomer if more than Xbp from end and less than Xbp from begin
    seqLen = len(sequence)
    dist2beg = candidate.beg
    dist2end = seqLen - candidate.end

    ## Filter out monomer if more than Xbp from end and less than Xbp from begin
    seqLen = len(sequence)
    dist2beg = candidate.beg
    dist2end = seqLen - candidate.end # Threshold at 50 as the TSD is on the 5' of the INS and TSD typically do not have sizes over this value

    if (dist2end <= 50) and (dist2beg >= 15):
        polyA = True

    else:
        polyA = False 

    return polyA, candidate

def search4polyT(sequence):
    '''
    Search for poly(T) at sequence begin 
    '''
    ## Configuration for monomere search:
    windowSize = 10
    maxWindowDist = 1
    minMonomerSize = 10
    minPurity = 80

    ## Seach poly(T) monomers
    targetMonomer = 'T'
    monomersT = sequences.find_monomers(sequence, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

    if not monomersT:
        return False, None

    ## Select monomer closest to sequence beg        
    candidate = monomersT[0]

    ## Filter out monomer if more than Xbp from end and less than Xbp from begin
    seqLen = len(sequence)
    dist2beg = candidate.beg  # Threshold at 50 as the TSD is on the 5' of the INS and TSD typically do not have sizes over this value
    dist2end = seqLen - candidate.end

    if (dist2beg <= 50) and (dist2end >= 15):
        polyT = True

    else:
        polyT = False 

    return polyT, candidate

def annotate_tsd(VCF, reference, retro_annot):
    '''
    '''
    ## Load reference using pysam
    reference = pysam.FastaFile(reference)

    ## For each variant
    for variant in VCF.variants:

        insId = variant.ID + '_' + variant.chrom + ':' + str(variant.pos) 

        ## Check for TSD downstream insertion position
        ## -------[TSD]INS[TARGET]-------
        targetCoord = variant.chrom + ':' + str(variant.pos) + '-' + str(variant.pos + 40)
        targetSeq = reference.fetch(region=targetCoord)
        tsdSeq_down, tsdLen_down, retro_annot[insId]['INS_TRIMMED'] = search4tsd(targetSeq, retro_annot[insId]['INS_TRIMMED'])

        ## Check for TSD upstream insertion position
        ## -------[TARGET]INS[TSD]-------
        targetCoord = variant.chrom + ':' + str(variant.pos - 40) + '-' + str(variant.pos)
        targetSeq = reference.fetch(region=targetCoord)
        tsdSeq_up, tsdLen_up, trimmed = search4tsd(targetSeq[::-1], retro_annot[insId]['INS_TRIMMED'][::-1])
        retro_annot[insId]['INS_TRIMMED'] = trimmed[::-1]

        ## Total TSD
        tsdSeq = tsdSeq_up[::-1] + tsdSeq_down
        tsdLen = tsdLen_up + tsdLen_down
        
        if tsdLen > 0:
            retro_annot[insId]['TSD_SEQ'] = tsdSeq
            retro_annot[insId]['TSD_LEN'] = tsdLen

def search4tsd(targetSeq, insert):
    '''
    '''
    ## Set to upper case to avoid format inconsistencies
    targetSeq = targetSeq.upper()
    insert = insert.upper()

    ## Set up parameters
    nbMM = 0
    maxMM = 2
    tsdLen = 0

    ## Pairwise compare
    for i in range(0, len(targetSeq)):
        if i >= len(insert):
            break 

        if targetSeq[i] != insert[i]:
            nbMM += 1
        
        if nbMM > maxMM:
            break

        tsdLen += 1
    
    ## Discard TSDs equal to the max number of allowed mismatches
    if tsdLen <= (maxMM + 1):
        tsdLen = 0
        tsdSeq = ''
        insert = insert
    else: 
        tsdSeq  = insert[0:tsdLen]
        insert = insert[tsdLen:]

    return tsdSeq, tsdLen, insert

def annotate_endomotif(VCF, reference, retro_annot):
    '''
    TODO
    '''
    ## Load reference using pysam
    reference = pysam.FastaFile(reference)

    ## For each variant
    for variant in VCF.variants:

        insId = variant.ID + '_' + variant.chrom + ':' + str(variant.pos) 

        ## Check for TSD downstream insertion position
        ## -------[TSD]INS[TARGET]-------
        targetCoord = variant.chrom + ':' + str(variant.pos) + '-' + str(variant.pos + 40)
        targetSeq = reference.fetch(region=targetCoord)

        print('ENDOMOTIF? ', insId, '+', retro_annot[insId]['STRAND'], targetSeq)


        ## Check for TSD upstream insertion position
        ## -------[TARGET]INS[TSD]-------
        targetCoord = variant.chrom + ':' + str(variant.pos - 40) + '-' + str(variant.pos)
        targetSeq = reference.fetch(region=targetCoord)
        
        #tsdSeq_up, tsdLen_up, trimmed = search4tsd(targetSeq[::-1], retro_annot[insId]['INS_TRIMMED'][::-1])
        
        print('ENDOMOTIF? ', insId, '-', retro_annot[insId]['STRAND'], targetSeq)


def rev2forward(retro_annot):
    '''
    '''
    for insId, annot in retro_annot.items():
        if annot['STRAND'] == '-':
            retro_annot[insId]['INS_TRIMMED'] = sequences.rev_complement(annot['INS_TRIMMED'])

def processed_pseudogene_annot(retro_annot, PAF_ref, exonsBinDb):
    '''
    '''
    ## 1. Filter out hits with MAPQ == 0
    PAF_ref.alignments = [hit for hit in PAF_ref.alignments if hit.MAPQ > 0]

    ## 2. Skip no hits on the reference
    if not PAF_ref.alignments:
        return retro_annot

    ## 3. Create hit chain
    chain = PAF_ref.chain(10, 25)

    ## 4. Intersects hits with annotated exons database
    for hit in chain.alignments:

        hit.geneNames = []

        ## Discard hits over references not included in the bin database
        if hit.tName not in list(exonsBinDb.keys()):
            continue

        ## Do intersection
        overlaps = exonsBinDb[hit.tName].collect_interval(hit.tBeg, hit.tEnd, ['EXON'])

        ## Collect gene names
        hit.geneNames = list(set([i[0].optional['geneName'] for i in overlaps]))
        
    ## 5. Filter out insertions without hits in exons
    nbHitsInExons = len([hit.geneNames for hit in chain.alignments if hit.geneNames])

    if nbHitsInExons == 0:
        return retro_annot
        
    ## 6. Define parent gene aggregating overlaps from all hits
    geneNames = [hit.geneNames for hit in chain.alignments if hit.geneNames]
    geneNames = [i for subList in geneNames for i in subList] ## Flatten the geneName list
    geneNamesCounts = dict(Counter(geneNames))

    ## 7. Select the gene name with the maximum number of assigned exons
    geneName = max(geneNamesCounts, key=geneNamesCounts.get)
    nbHitsInExons = geneNamesCounts[geneName]

    ## 8. Annotate insertion structure 
    retro_annot['CONFORMATION'] = annotate_structure_pseudogene(PAF_ref, retro_annot['STRAND'],retro_annot['INS_TRIMMED'])

    ## 9. Annotate PSD poly(A)
    targetMonomer = 'A'
    windowSize = 10
    maxWindowDist = 1
    minMonomerSize = 10
    minPurity = 80

    monomers = sequences.find_monomers_backward(retro_annot['INS_TRIMMED'], targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

    if monomers:
        retro_annot['POLYA_SEQ'] = monomers[-1].seq
        retro_annot['POLYA_LEN'] = len(retro_annot['POLYA_SEQ'])
        retro_annot['INS_TRIMMED'] = retro_annot['INS_TRIMMED'][:monomers[-1].beg]

    ## 10. Save PSD annotation info
    retro_annot['DTYPE_N'] = 'PSD'
    retro_annot['SRC_GENE'] = geneName
    retro_annot['NB_EXONS'] = nbHitsInExons
    retro_annot['PSD_LEN'] = len(retro_annot['INS_TRIMMED'])


    return retro_annot


def annotate_structure_pseudogene(PAF_ref, strand, seq):
    '''
    '''
    configuration = ''
    prev_operation = ''

    for hit in PAF_ref.alignments:

        ## Define operation
        if (strand == '-') and (hit.strand == '+'):
            operation = 'REV'
        elif (strand == '-') and (hit.strand == '-'):
            operation = 'FOR'
        elif (strand == '+') and (hit.strand == '+'):
            operation = 'FOR'
        elif (strand == '+') and (hit.strand == '-'):
            operation = 'REV'

        ## Skip if same operation as previous
        if operation == prev_operation:
            continue

        ## Add operation to the configuration
        if configuration == '':
            configuration = operation
        else:
            configuration = configuration + '+' + operation

        ## update previous operation
        prev_operation = operation
    
    configuration = configuration + '+POLYA'

    return configuration

def search4transductions_3prime(retro_annot, exonsBinDb, reference, outDir):
    '''
    '''
    ## For each insertion
    for insId, annot in retro_annot.items():

        ## Filter out insertions with already assigned insert type
        if annot['DTYPE_N'] != 'UNK':
            continue

        ## Annotate and trim transduced sequence plus polyAs
        retro_annot[insId] = annotate_3prime_td(annot['INS_TRIMMED'], annot)

    ## Align transduced sequences against the reference or database of sequences downstream of potential source elements
    ## Note: database sample-specific for cancer genomes vs catalogue for germline polymorphisms
    retro_annot = align_3prime_td(retro_annot, exonsBinDb, reference, outDir)

def annotate_3prime_td(seq, annot):
    '''
    TO DO Tomorrow, Align transduced sequences on the reference and do assignations to source regions
    '''
    ## Call monomers
    targetMonomer = 'A'
    windowSize = 10
    maxWindowDist = 1
    minMonomerSize = 10
    minPurity = 80

    monomers = sequences.find_monomers_backward(seq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

    ## Filter internal monomers
    maxDist2end = 20
    minInternalMonomerSize = 15

    monomers = sequences.filter_no3prime_monomers(monomers, seq, maxDist2end, minInternalMonomerSize)

    ## A) No monomers. No polyA 
    if len(monomers) == 0:
        annot['INS_TRIMMED'] = seq

    ## B) One monomer. Single polyA
    elif len(monomers) == 1:
        annot['POLYA_LEN'] = str(monomers[0].length())
        annot['INS_TRIMMED'] = seq[0:monomers[0].beg]
        annot['POLYA_SEQ'] = seq[monomers[0].beg:]

    ## C) Two or more monomers. Partnered/multi-partnered transduction
    else:
        annot['DTYPE_N'] = 'partnered'
        annot['POLYA_LEN'] = ','.join([str(monomer.length()) for monomer in monomers])
        annot['3PRIME_TD_LEN'] = ','.join([str(len(seq[monomerA.end:monomerB.beg])) for monomerA, monomerB in zip(monomers, monomers[1:])])
        annot['INS_TRIMMED'] = seq[0:monomers[0].beg]
        annot['POLYA_SEQ'] = ','.join([monomer.seq for monomer in monomers])
        annot['3PRIME_TD_SEQ'] = ','.join([seq[monomerA.end:monomerB.beg] for monomerA, monomerB in zip(monomers, monomers[1:])])
        annot['3PRIME_NB_TD'] = len(monomers) - 1 

    return annot

def align_3prime_td(retro_annot, exonsBinDb, reference, outDir):
    '''
    TO DO: MERGE transduced hits from same source region

    NOTE: Not all information included in the dictionary
    '''
    ## Create fasta with 3prime transduction sequences
    FASTA_path = create_fasta_3prime_td(retro_annot, outDir)

    ## Align transduced sequences into the reference
    #PAF_path = '/Users/brodriguez/Research/Projects/1KGP_LR/SV_annot/resubmission/annot_INS/tmp/td3prime2ref.paf'
    SAM_path = alignment.alignment_bwa(FASTA_path, reference, 'td3prime2ref', 2, outDir)     
    PAF_path = alignment.sam2paf(SAM_path, 'td3prime2ref', outDir)
        
    ## Read alignments and create chain
    PAF_ref = formats.PAF()
    PAF_ref.read(PAF_path)

    ## Filter out hits with MAPQ == 0
    PAF_ref.alignments = [hit for hit in PAF_ref.alignments if hit.MAPQ > 0]
    
    ## Discard if no hit passes the filters
    if not PAF_ref.alignments:
        return
    
    ## Group alignments by transduction id
    refHitsTd = group_alignments(PAF_ref)

    ## For each insertion, add transduction alignment coordinates
    for insId in retro_annot:
        if ('3PRIME_NB_TD' not in retro_annot[insId]) or (retro_annot[insId]['3PRIME_NB_TD'] == 0):
            continue

        tdCoords = []
        chain = None

        for i in range(1, (retro_annot[insId]['3PRIME_NB_TD'] + 1)):
            tdId = insId + '_' + str(i)

            if tdId not in refHitsTd:
                tdCoords.append('UNK')
                continue

            ## create chain
            chain = refHitsTd[tdId].chain(10, 25)

            ## Transduced sequence able to trace back to source locus
            if chain.perc_query_covered() > 80:
                ref = chain.alignments[0].tName
                beg, end = chain.interval_template()
                coord = ref + ':' + str(beg) + '-' + str(end)
                tdCoords.append(coord)

                ## Check if transduced sequence encompass exonic sequences
                minPercOverlap = 80

                for refHit in chain.alignments:

                    refHit.geneNames = []
                    
                    if refHit.tName not in list(exonsBinDb.keys()):
                        continue

                    exonOverlaps = exonsBinDb[refHit.tName].collect_interval(refHit.tBeg, refHit.tEnd, ['EXON'])

                    ## Filter overlaps
                    filtered_exonOverlaps = [i for i in exonOverlaps if i[2] >= minPercOverlap]

                    ## Collect gene names
                    if filtered_exonOverlaps:
                        refHit.geneNames = list(set([i[0].optional['geneName'] for i in filtered_exonOverlaps]))

                ## Collect geneName info
                geneNames = [hit.geneNames for hit in chain.alignments if hit.geneNames]
                geneNames = list(set([i for subList in geneNames for i in subList])) ## Flatten the geneName list
                geneNames = ','.join(geneNames)
            
                if geneNames:
                    retro_annot[insId]['3PRIME_MOB_GENE'] = geneNames

            else:
                ## Candidate transduced sequence does not align on the reference. 
                tdCoords.append('UNK')

        ## Convert into string and add to the annot dictionary
        retro_annot[insId]['3PRIME_TD_COORD'] = ','.join(tdCoords)

        if chain is not None:
            retro_annot[insId]['3PRIME_TD_MAPQ'] = ','.join([str(hit.MAPQ) for hit in chain.alignments])

    return retro_annot
    ## Cleanup
    #unix.rm([FASTA_path, SAM_path, PAF_path, outDir + '/align.err'])

def create_fasta_3prime_td(retro_annot, outDir):
    '''
    '''
    ## Initialize fasta object
    fasta = formats.FASTA()

    ## Add td sequences to fasta object
    for insId, annot in retro_annot.items():
        if ('3PRIME_NB_TD' not in annot) or (annot['3PRIME_NB_TD'] == 0):
            continue

        tdSeqs = annot['3PRIME_TD_SEQ'].split(',')

        for i in range(annot['3PRIME_NB_TD']):

            tdId = insId + '_' + str(i + 1)
            td = tdSeqs[i]
            fasta.seqDict[tdId] = td

    ## Write sequences to fasta file
    FASTA_path = outDir + '/3prime_transductions.fa'
    fasta.write(FASTA_path)

    return FASTA_path

def create_fasta_nested_ins(retro_annot, outDir):
    '''
    '''
    ## Initialize fasta object
    fasta = formats.FASTA()

    ## Add td sequences to fasta object
    for insId, annot in retro_annot.items():

        if ('3PRIME_TD_COORD' in annot) and (annot['3PRIME_TD_COORD'] == 'UNK'):
            fasta.seqDict[insId] = annot['3PRIME_TD_SEQ']

    ## Write sequences to fasta file
    FASTA_path = outDir + '/nested_ins_candidates.fa'
    fasta.write(FASTA_path)

    return FASTA_path

def annot_retroelement(retro_annot, exonsBinDb, repeatsDb, TRF, reference, outDir):
    '''
    '''
    ## 1. Align candidate retroelement inserts against retroelements database 
    PAFs_repeats = alignTrimmedIns2repeats(retro_annot, repeatsDb, 0, 15, 'trimmed2repeats', outDir)

    ## 2. Align candidate retrotransposition inserts against reference genome
    PAFs_ref = alignTrimmedIns2ref(retro_annot, reference, 10, 15, 'trimmed2ref', outDir)

    ## 3. Filter reference hits overlapping with repeat hits
    PAFs_ref_pass = filter_overlapping_hits(PAFs_ref, PAFs_repeats, 10)

    ## 4. Merge retro and reference hits
    PAFs_merged = merge_hits(PAFs_repeats, PAFs_ref_pass)

    ## 5. Create hit chains, one per insertion
    for insId in PAFs_repeats.keys():
        retro_annot[insId]['CHAIN_repeats'] = PAFs_repeats[insId].chain(300, 25) ## NOTE: parameters to tune, specially the distance threshold

    for insId in PAFs_merged.keys():
        retro_annot[insId]['CHAIN_merged'] = PAFs_merged[insId].chain(300, 25) ## NOTE: parameters to tune, specially the distance threshold

    ## 6. Orphan transduction annotation  
    for insId in PAFs_merged.keys():
        retro_annot[insId] = orphan_td_annot(retro_annot[insId], retro_annot[insId]['CHAIN_merged'])

    ## 7. Templated insertions
    for insId in PAFs_merged.keys():
        retro_annot[insId] = templated_ins(retro_annot[insId], retro_annot[insId]['CHAIN_merged'])

    ## 8. 5' transductions.
    for insId in PAFs_repeats.keys():

        PAF_ref = PAFs_ref[insId] if insId in PAFs_ref else formats.PAF()
        retro_annot[insId] = transductions_5prime(retro_annot[insId], retro_annot[insId]['CHAIN_repeats'], PAF_ref, exonsBinDb)

    ## 9. 3' transduction refinement.
    for insId in PAFs_repeats.keys():
        PAF_ref = PAFs_ref[insId] if insId in PAFs_ref else formats.PAF()
        transductions_3prime_hits(retro_annot[insId], retro_annot[insId]['CHAIN_repeats'], PAF_ref, exonsBinDb)

    ## 10. 5' transductions refinement
    retro_annot = transductions_5prime_refinement(retro_annot, reference, exonsBinDb, 10, 15, outDir)
    
    ## 11. Reclassify partnered with the transduced sequencing aligning at the insertion position. Probably are insertions containing templated insertions
    retro_annot = reclassify_5prime_partnered_ipos(retro_annot)

    ## 12. Family assignation
    for insId in PAFs_repeats.keys():
        retro_annot[insId] = family_annot_insert(retro_annot[insId], retro_annot[insId]['CHAIN_repeats'])

    ## 13. Solo retrotransposon annotation
    for insId in PAFs_repeats.keys():
        retro_annot[insId] = solo_annot(retro_annot[insId], retro_annot[insId]['CHAIN_repeats'])

    ## 14. Nested insertion annotation
    retro_annot = nested_annot(retro_annot, repeatsDb, 0, 15, 'td2repeats', outDir)

    ## 15. Compute % resolved for every insertion class
    for insId in retro_annot:
        retro_annot[insId] = compute_perc_resolved(retro_annot[insId])

    ## 16. Reclasify as solo those transductions not traced to their source element
    for insId in retro_annot:
        retro_annot[insId] = reclassify_partnered_not_traced(retro_annot[insId])

    ## 17. Filter out insertions whose type could not be determined
    retro_annot_pass = dict([(insId, retro_annot[insId]) for insId in retro_annot if (retro_annot[insId]['DTYPE_N'] != 'UNK')])
    retro_annot_fail = dict([(insId, retro_annot[insId]) for insId in retro_annot if (retro_annot[insId]['DTYPE_N'] == 'UNK')])
    print('NB_RETRO_FAIL_ITYPE: ', len(retro_annot_fail))

    ## 18. Filter out based on % resolved 
    ## NOTE: SVAs still not perfectly annotated so lower % resolved in comparison to L1/Alus
    retro_annot_pass, filter_fail = filter_perc_resolved(retro_annot_pass)
    retro_annot_fail.update(filter_fail)
    print('NB_RETRO_FAIL_PERC_RESOLVED: ', len(retro_annot_fail))

    ## 19. Filter out partnered with the transduced sequencing aligning at the insertion position. Probably are tandem DUPs
    retro_annot_pass, filter_fail = filter_partnered_ipos(retro_annot_pass)
    retro_annot_fail.update(filter_fail)
    print('NB_RETRO_FAIL_PARTNERED_IPOS: ', len(retro_annot_fail))

    ## 20. Filter out non-sva insertions with hexamer annotation
    retro_annot_pass, filter_fail = filter_no_sva_hexamer(retro_annot_pass)
    retro_annot_fail.update(filter_fail)
    print('NB_RETRO_FAIL_HEXAMER: ', len(retro_annot_fail))

    ## 21. Filter out non-sva insertions with hexamer annotation
    retro_annot_pass, filter_fail = filter_PSD_MAST2(retro_annot_pass)
    retro_annot_fail.update(filter_fail)
    print('NB_RETRO_FAIL_PSD_MAST2: ', len(retro_annot_fail))

    ## 22. Filter out PSD insertions with templated insertions
    retro_annot_pass, filter_fail = filter_PSD_TEMPLATED(retro_annot_pass)
    retro_annot_fail.update(filter_fail)
    print('NB_RETRO_FAIL_PSD_TEMPLATED: ', len(retro_annot_fail))

    ## 23. Filter out ERVK insertions as they should not be detected in this module
    retro_annot_pass, filter_fail = filter_HERVK(retro_annot_pass)
    retro_annot_fail.update(filter_fail)
    print('NB_RETRO_FAIL_HERVK: ', len(retro_annot_fail))

    ## 24. Annotate insertion conformation and structure for those passing the filters
    for insId in retro_annot_pass:

        ## Initialize insertion status as not canonical
        retro_annot_pass[insId]['NOT_CANONICAL'] = True

        ## A) L1 and Alu
        if ('FAM_N' in retro_annot_pass[insId]) and (retro_annot_pass[insId]['FAM_N'] in ['L1', 'Alu']):
            retro_annot_pass[insId] = annotate_structure(retro_annot_pass[insId])

        ## B) SVA
        elif ('FAM_N' in retro_annot_pass[insId]) and (retro_annot_pass[insId]['FAM_N'] == 'SVA'):
            ## Structure
            retro_annot_pass[insId] = annotate_structure_sva(retro_annot_pass[insId])

            ## VNTR
            retro_annot_pass[insId] = annotate_VNTR(retro_annot_pass[insId])
        
        ## Determine if not canonical conformation
        if ('FAM_N' in retro_annot_pass[insId]) and ('CONFORMATION' in retro_annot_pass[insId]):

            ## L1
            if (retro_annot_pass[insId]['FAM_N'] == 'L1'):
                canonicalList = ['FOR+POLYA', 'TRUN+FOR+POLYA', 'TRUN+REV+DEL+FOR+POLYA', 'TRUN+REV+DUP+FOR+POLYA', 'TRUN+REV+BLUNT+FOR+POLYA', 'FOR+POLYA+TD+POLYA', 'TRUN+FOR+POLYA+TD+POLYA', 'TRUN+REV+DEL+FOR+POLYA+TD+POLYA', 'TRUN+REV+DUP+FOR+POLYA+TD+POLYA', 'TRUN+REV+BLUNT+FOR+POLYA+TD+POLYA', 'TD+FOR+POLYA', 'REV+DEL+FOR+POLYA']

            ## Alu
            elif (retro_annot_pass[insId]['FAM_N'] == 'Alu'):
                canonicalList = ['TRUN+FOR+POLYA', 'FOR+POLYA']

            ## SVA
            else:
                canonicalList = ['Hexamer+Alu-like+VNTR+SINE-R+POLYA', 'Alu-like+VNTR+SINE-R+POLYA', 'VNTR+SINE-R+POLYA', 'SINE-R+POLYA', 'Hexamer+Alu-like+VNTR+SINE-R+POLYA+TD+POLYA', 'Alu-like+VNTR+SINE-R+POLYA+TD+POLYA', 'VNTR+SINE-R+POLYA+TD+POLYA', 'SINE-R+POLYA+TD+POLYA', 'TD+Hexamer+Alu-like+VNTR+SINE-R+POLYA', 'MAST2+VNTR+SINE-R+POLYA', 'MAST2+VNTR+SINE-R+POLYA+TD+POLYA', 'TD+MAST2+VNTR+SINE-R+POLYA']

            CANONICAL = is_canonical(retro_annot_pass[insId]['CONFORMATION'], canonicalList)

            if CANONICAL:
                retro_annot_pass[insId]['NOT_CANONICAL'] = False

    return retro_annot_pass, retro_annot_fail

def reclassify_partnered_not_traced(annot):
    '''
    '''
    if '3PRIME_TD_COORD' not in annot:
        return annot
    
    sourceCoord = annot['3PRIME_TD_COORD'].split(',')[0]

    ## Reclassify as solo and remove transduction information
    if sourceCoord == 'UNK':

        ## Redefine itype
        annot['DTYPE_N'] = 'solo'

        ## Report last polyA tail
        annot['POLYA_LEN'] = annot['POLYA_LEN'].split(',')[-1]
        annot['POLYA_SEQ'] = annot['POLYA_SEQ'].split(',')[-1]

        ## Remove TD related info from the conformation
        if 'CONFORMATION' in annot:
            annot['CONFORMATION'] = annot['CONFORMATION'].split('+TD+POLYA')[0]

        ## Delete transduction related fields
        del annot['3PRIME_NB_TD']
        del annot['3PRIME_TD_LEN']
        del annot['3PRIME_TD_COORD']
        del annot['3PRIME_TD_SEQ']

    return annot

def filter_no_sva_hexamer(retro_annot):
    '''
    '''
    filter_pass = {} 
    filter_fail = {}

    for insId in retro_annot:

        if 'HEXAMER_LEN' not in retro_annot[insId]:
            filter_pass[insId] = retro_annot[insId]
            continue

        if ('FAM_N' in retro_annot[insId]) and (retro_annot[insId]['FAM_N'] == 'SVA'):
            filter_pass[insId] = retro_annot[insId]
            continue

        elif (retro_annot[insId]['DTYPE_N'] == 'PSD') or (retro_annot[insId]['DTYPE_N'] == 'orphan'):
            filter_pass[insId] = retro_annot[insId]
            del filter_pass[insId]['HEXAMER_LEN']
            del filter_pass[insId]['HEXAMER_SEQ']
            continue

        else:
            filter_fail[insId] = retro_annot[insId]

    return filter_pass, filter_fail

def filter_PSD_MAST2(retro_annot):
    '''
    '''
    filter_pass = {} 
    filter_fail = {}

    for insId in retro_annot:

        if (retro_annot[insId]['DTYPE_N'] == 'PSD') and (retro_annot[insId]['SRC_GENE'] == 'MAST2'):
            filter_fail[insId] = retro_annot[insId]
        else:
            filter_pass[insId] = retro_annot[insId]

    return filter_pass, filter_fail


def filter_PSD_TEMPLATED(retro_annot):
    '''
    '''
    filter_pass = {} 
    filter_fail = {}

    for insId in retro_annot:

        if (retro_annot[insId]['DTYPE_N'] == 'PSD') and (('5PRIME_TEMP_LEN' in retro_annot[insId]) or ('3PRIME_TEMP_LEN' in retro_annot[insId])):
            filter_fail[insId] = retro_annot[insId]
        else:
            filter_pass[insId] = retro_annot[insId]

    return filter_pass, filter_fail

def filter_HERVK(retro_annot):
    '''
    '''
    filter_pass = {} 
    filter_fail = {}

    for insId in retro_annot:

        if ('FAM_N' in retro_annot[insId]) and (retro_annot[insId]['FAM_N'] == 'HERVK'):
            filter_fail[insId] = retro_annot[insId]
        else:
            filter_pass[insId] = retro_annot[insId]

    return filter_pass, filter_fail


def is_canonical(conformation, canonicalList):
    '''
    '''
    CANONICAL = False
    canonicalDict = {}

    for conf in canonicalList:
        canonicalDict[conf] = 1

    if conformation in canonicalDict:
        CANONICAL = True
    
    return CANONICAL

def annotate_hexamer(VCF, retro_annot, TRF_out):
    '''
    '''
    ### Load TRF calls
    TRF = formats.TRF()
    TRF.read(TRF_out)

    ## For each variant
    for variant in VCF.variants:

        insId = variant.ID + '_' + variant.chrom + ':' + str(variant.pos) 

        ## Skip if not retrotransposon candidate or not TRF call
        if (insId not in retro_annot) or (insId not in TRF.callsDict):
            retro_annot[insId]['INS_TRIMMED'] = variant.ref
            continue
        
        hexamerSeq, retro_annot[insId]['INS_TRIMMED'] = TRF.callsDict[insId].SVA_hexamer(variant.ref, retro_annot[insId]['STRAND'], 10, 80)

        ## Hexamer found
        if hexamerSeq is not None:
            retro_annot[insId]['HEXAMER_SEQ'] = hexamerSeq
            retro_annot[insId]['HEXAMER_LEN'] = len(hexamerSeq) 

def filter_perc_resolved(retro_annot):
    '''
    '''
    filter_pass = {} 
    filter_fail = {}

    for insId in retro_annot:

        ## A) Orphan or PSD
        if (retro_annot[insId]['DTYPE_N'] == 'orphan') or (retro_annot[insId]['DTYPE_N'] == 'PSD'):

            if retro_annot[insId]['PERC_RESOLVED'] >= 75:
                filter_pass[insId] = retro_annot[insId]        

            else:
                filter_fail[insId] = retro_annot[insId] 
            
        ## B) Chimera
        elif (retro_annot[insId]['DTYPE_N'] == 'chimera'):
            filter_pass[insId] = retro_annot[insId]        

        ## C) solo or partnered without assigned family
        elif ((retro_annot[insId]['DTYPE_N'] == 'solo') or (retro_annot[insId]['DTYPE_N'] == 'partnered')) and ('FAM_N' not in retro_annot[insId]):

            filter_fail[insId] = retro_annot[insId]                

        ## D) solo or partnered SVA
        elif ((retro_annot[insId]['DTYPE_N'] == 'solo') or (retro_annot[insId]['DTYPE_N'] == 'partnered')) and (retro_annot[insId]['FAM_N'] == 'SVA'):

            if retro_annot[insId]['PERC_RESOLVED'] >= 50:
                filter_pass[insId] = retro_annot[insId]        

            else:
                filter_fail[insId] = retro_annot[insId] 

        ## E) solo or partnered other element (Alu/L1)
        elif ((retro_annot[insId]['DTYPE_N'] == 'solo') or (retro_annot[insId]['DTYPE_N'] == 'partnered')):

            if retro_annot[insId]['PERC_RESOLVED'] >= 75:
                filter_pass[insId] = retro_annot[insId]        

            else:
                filter_fail[insId] = retro_annot[insId] 

    return filter_pass, filter_fail

def filter_partnered_ipos(retro_annot):
    '''
    '''
    filter_pass = {} 
    filter_fail = {}

    for insId in retro_annot:

        ## Pass if not partnered transduction
        if retro_annot[insId]['DTYPE_N'] != 'partnered':
            filter_pass[insId] = retro_annot[insId]
            continue

        ## Collect coordinates for insertion interval 
        tmp = insId.split('_')[1]
        iChr, iPos = tmp.split(':')
        iBeg = int(iPos) - 1000
        iEnd = int(iPos) + 1000

        ## Collect coordinates for transduced regions
        tdCoords = ''

        if '3PRIME_TD_COORD' in retro_annot[insId]: 
            tdCoords = retro_annot[insId]['3PRIME_TD_COORD'] + ','

        if '5PRIME_TD_COORD' in retro_annot[insId]: 
            tdCoords = tdCoords + retro_annot[insId]['5PRIME_TD_COORD']

        ## Assess if transduced region coordinates overlap with the insertion interval
        overlap = False

        for tdCoord in tdCoords.split(','):
            if (tdCoord == 'UNK') or (tdCoord == ''):
                continue
            
            tdChr, tdInterval = tdCoord.split(':')
            tdBeg, tdEnd = tdInterval.split('-')

            if (iChr == tdChr) and (gRanges.overlap(iBeg, iEnd, int(tdBeg), int(tdEnd))[0]):
                overlap = True

        ## Filter out partnered transduction if overlap found (likely tandem duplication)
        if overlap:
            filter_fail[insId] = {}
            
        else:
            filter_pass[insId] = retro_annot[insId]

    return filter_pass, filter_fail

def reclassify_5prime_partnered_ipos(retro_annot):
    '''
    '''
    for insId in retro_annot:

        ## Pass if not 5' partnered transduction
        if (retro_annot[insId]['DTYPE_N'] != 'partnered') or ('5PRIME_TD_COORD' not in retro_annot[insId]):
            continue

        if (retro_annot[insId]['5PRIME_NB_TD'] > 1):
            continue

        ## Collect coordinates for insertion interval 
        tmp = insId.split('_')[1]
        iChr, iPos = tmp.split(':')
        iBeg = int(iPos) - 1000
        iEnd = int(iPos) + 1000

        ## Collect coordinates for transduced regions
        tdCoords = retro_annot[insId]['5PRIME_TD_COORD']

        ## Assess if transduced region coordinates overlap with the insertion interval
        overlap = False

        for tdCoord in tdCoords.split(','):
            if (tdCoord == 'UNK') or (tdCoord == ''):
                continue
            
            tdChr, tdInterval = tdCoord.split(':')
            tdBeg, tdEnd = tdInterval.split('-')

            if (iChr == tdChr) and (gRanges.overlap(iBeg, iEnd, int(tdBeg), int(tdEnd))[0]):
                overlap = True

        ## Reclassify partnered transductions as potential templated insertions
        if overlap:

            if ('3PRIME_TD_COORD' not in retro_annot[insId]):
                retro_annot[insId]['DTYPE_N'] = 'UNK'
            
            retro_annot[insId]['5PRIME_TEMP_SEQ'] = retro_annot[insId]['5PRIME_TD_SEQ']
            retro_annot[insId]['5PRIME_TEMP_LEN'] = int(retro_annot[insId]['5PRIME_TD_LEN'])
            retro_annot[insId]['5PRIME_TEMP_COORD'] = retro_annot[insId]['5PRIME_TD_COORD']

            del retro_annot[insId]['5PRIME_TD_SEQ']
            del retro_annot[insId]['5PRIME_TD_LEN']
            del retro_annot[insId]['5PRIME_TD_COORD']
            del retro_annot[insId]['5PRIME_TD_MAPQ']
            del retro_annot[insId]['5PRIME_NB_TD']

    return retro_annot

def transductions_5prime_refinement(retro_annot, reference, exonsBinDb, minMAPQ, minHitLen, outDir):
    '''
    '''

    ## 1. Collect candidate 5' transduction sequences
    FASTA_5prime = formats.FASTA() 

    for insId in retro_annot:

        ## Skip insertion without a chain or trimmed inserted sequence available
        if ('CHAIN_merged' not in retro_annot[insId]) or ('INS_TRIMMED' not in retro_annot[insId]):
            continue

        ## Skip insertions with already 5' transduction calls OR orphan OR pseudogene
        if ('5PRIME_TD_LEN' in retro_annot[insId]) or (retro_annot[insId]['DTYPE_N'] == 'orphan') or (retro_annot[insId]['DTYPE_N'] == 'PSD'):
            continue

        ## Skip insertion if not long-enough sequence on the 5' to realign
        beg = retro_annot[insId]['CHAIN_merged'].interval()[0]

        if beg < 15:
            continue

        ## Add to the fasta sequence to realign
        seq2realign = retro_annot[insId]['INS_TRIMMED'][0:beg]
        FASTA_5prime.seqDict[insId] = seq2realign

    ## 2. Write candidate transduction sequences into a fasta
    FASTA_path = outDir + '/5prime_td_candidates.fa'        
    FASTA_5prime.write(FASTA_path)

    ## 3. Align sequences with bwa mem
    #PAF_path = '/Users/brodriguez/Research/Projects/1KGP_LR/SV_annot/resubmission/annot_INS/tmp/5prime_td_candidates.paf'
    SAM_path = alignment.alignment_bwa(FASTA_path, reference, '5prime_td_candidates', 2, outDir)
    PAF_path = alignment.sam2paf(SAM_path, '5prime_td_candidates', outDir)
        
    ## 4. Read alignments file
    PAF_hits = formats.PAF()
    PAF_hits.read(PAF_path)

    ## 3. Filter hits based on MAPQ
    PAF_hits.alignments = [hit for hit in PAF_hits.alignments if hit.MAPQ >= minMAPQ]

    ## 4. Filter hits shorter than Xbp
    PAF_hits.alignments = [hit for hit in PAF_hits.alignments if hit.alignmentLen() >= minHitLen]

    ## 5. Generate single PAF objects per inserted sequence:
    PAFs_repeats = group_alignments(PAF_hits)

    ## 5. Generate single PAF objects per inserted sequence:
    PAFs_hits = group_alignments(PAF_hits)

    ## 6. For each insertion, add transduction alignment coordinates
    for insId in PAFs_hits:

        ## create chain
        chain = PAFs_hits[insId].chain(10, 25)
        
        ## Transduced sequence able to trace back to source locus
        if chain.perc_query_covered() > 80:

            ## Transduction coordinates and MAPQ
            retro_annot[insId]['DTYPE_N'] = 'partnered'
            retro_annot[insId]['5PRIME_TD_LEN'] = ','.join([str(hit.qEnd - hit.qBeg) for hit in chain.alignments])
            retro_annot[insId]['5PRIME_TD_SEQ'] = ','.join([retro_annot[insId]['INS_TRIMMED'][hit.qBeg:hit.qEnd] for hit in chain.alignments])
            retro_annot[insId]['5PRIME_NB_TD'] = len(chain.alignments)
            retro_annot[insId]['5PRIME_TD_COORD'] = ','.join([hit.tName + ':' + str(hit.tBeg) + '-' + str(hit.tEnd) for hit in chain.alignments])
            retro_annot[insId]['5PRIME_TD_MAPQ'] = ','.join([str(hit.MAPQ) for hit in chain.alignments])

            ## Check if transduced sequence encompass exonic sequences
            minPercOverlap = 80

            for refHit in chain.alignments:

                refHit.geneNames = []

                if refHit.tName not in list(exonsBinDb.keys()):
                    continue

                exonOverlaps = exonsBinDb[refHit.tName].collect_interval(refHit.tBeg, refHit.tEnd, ['EXON'])                
    
                ## Filter overlaps
                filtered_exonOverlaps = [i for i in exonOverlaps if i[2] >= minPercOverlap]
                
                ## Collect gene names
                if filtered_exonOverlaps:
                    refHit.geneNames = list(set([i[0].optional['geneName'] for i in filtered_exonOverlaps]))

            ## Collect geneName info
            geneNames = [hit.geneNames for hit in chain.alignments if hit.geneNames]
            geneNames = list(set([i for subList in geneNames for i in subList])) ## Flatten the geneName list
            geneNames = ','.join(geneNames)
            
            if geneNames:
                retro_annot[insId]['5PRIME_MOB_GENE'] = geneNames
    
    return retro_annot

def annotate_structure_sva(annot):
    '''
    TO DO: I have noticed that the chain is incomplete. Somes hits are missing leading to lower % resolved rates. Can be related with consensus database or with the chaining process. Fix this issue first before implementing this function
    '''
    if not 'CHAIN_repeats' in annot:
        return annot
    
    ## Compute conformation avoiding redundancies
    conformationList = []

    for hit in annot['CHAIN_repeats'].alignments:
        feature = hit.tName.split('|')[-1]

        if feature not in conformationList:
            conformationList.append(feature)
            
    annot['CONFORMATION'] = '+'.join(conformationList)

    ## Compute extended conformation (TO DO)

    ## Add hexamer
    if 'HEXAMER_LEN' in annot:
        annot['CONFORMATION'] = 'Hexamer' + '+' + annot['CONFORMATION']

    ## Add 5' transduction
    if '5PRIME_NB_TD' in annot:
        annot['CONFORMATION'] = 'TD+' * int(annot['5PRIME_NB_TD']) + annot['CONFORMATION']

    ## Add 3' transduction
    if '3PRIME_NB_TD' in annot:
        annot['CONFORMATION'] = annot['CONFORMATION'] + int(annot['3PRIME_NB_TD']) * '+POLYA+TD'

    ## Add templated insertion
    if '5PRIME_TEMP_LEN' in annot:
        annot['CONFORMATION'] = 'TEMPLATED+' + annot['CONFORMATION']

    if '3PRIME_TEMP_LEN' in annot:
        annot['CONFORMATION'] = annot['CONFORMATION'] + '+TEMPLATED' 

    ## Add polyA
    if 'POLYA_LEN' in annot:
        annot['CONFORMATION'] = annot['CONFORMATION'] + '+POLYA'

    ## Add TSD
    #if 'TSD_LEN' in annot:
    #    annot['CONFORMATION'] = 'TSD+' + annot['CONFORMATION'] + '+TSD'

    return annot
    
def annotate_VNTR(annot):
    '''
    '''
    ## Collect all VNTRs
    VNTRs = []
    for hit in annot['CHAIN_repeats'].alignments:
        feature = hit.tName.split('|')[-1]
        if feature == 'VNTR':
            VNTRs.append(hit)

    if not VNTRs:
        return annot
    
    ## Determine VNTR length
    beg = VNTRs[0].qBeg
    end = VNTRs[-1].qEnd
    annot['VNTR_LEN'] = end - beg
    annot['VNTR_COORD'] = str(beg) + '-' + str(end)

    return annot


def annotate_structure(annot):
    '''
    '''
    if not 'CHAIN_repeats' in annot:
        return annot

    ## Annotate conformation for repetitive DNA element
    annot['CONFORMATION'], annot['CONFORMATION_EXT'] = annotate_conformation(annot['CHAIN_repeats'])

    ## Add 5' transduction
    if '5PRIME_NB_TD' in annot:
        annot['CONFORMATION'] = 'TD+' * int(annot['5PRIME_NB_TD']) + annot['CONFORMATION']

    ## Add 3' transduction
    if '3PRIME_NB_TD' in annot:
        annot['CONFORMATION'] = annot['CONFORMATION'] + int(annot['3PRIME_NB_TD']) * '+POLYA+TD'

    ## Add nested insertion
    if 'NESTED_FAM' in annot:
        annot['CONFORMATION'] = annot['CONFORMATION'] + '+POLYA+' + annot['NESTED_FAM']

    ## Add templated insertion
    if '5PRIME_TEMP_LEN' in annot:
        annot['CONFORMATION'] = 'TEMPLATED+' + annot['CONFORMATION']

    if '3PRIME_TEMP_LEN' in annot:
        annot['CONFORMATION'] = annot['CONFORMATION'] + '+TEMPLATED' 

    ## Add polyA
    if 'POLYA_LEN' in annot:
        annot['CONFORMATION'] = annot['CONFORMATION'] + '+POLYA'

    ## Add TSD
    #if 'TSD_LEN' in annot:
    #    annot['CONFORMATION'] = 'TSD+' + annot['CONFORMATION'] + '+TSD'

    return annot

def annotate_conformation(chain):
    '''
    '''
    ## Compute number of hits
    nbHits = len(chain.alignments)

    ## Check if there is 5' truncation
    minTrunLen = 10

    if int(chain.alignments[0].tBeg) > minTrunLen:
        trun5prime = 'TRUN'
        trun5prime_ext = 'TRUN' + '(' + str(chain.alignments[0].tBeg) + ',' + str(0) + '-' + str(chain.alignments[0].tBeg) + ')' 

    else:
        trun5prime = ''
        trun5prime_ext = ''

    ## Check if thre is 3' truncation
    minTrunLen = 10
    dist2end = int(chain.alignments[-1].tLen) - int(chain.alignments[-1].tEnd)

    if dist2end > minTrunLen:
        trun3prime = 'TRUN'
        trun3prime_ext = 'TRUN' + '(' + str(dist2end) + ',' + str(chain.alignments[-1].tEnd) + '-' + str(chain.alignments[-1].tLen) + ')' 

    else:
        trun3prime = ''
        trun3prime_ext = ''

    ## Single hit ##
    if nbHits == 1:
        conformation = 'FOR' if chain.alignments[0].strand == '+' else 'REV'
        conformation_ext = conformation + '(' + chain.alignments[0].alignmentString() + ')'    

        ## Add 5' truncation
        if trun5prime != '':
            conformation = trun5prime + '+' + conformation
            conformation_ext = trun5prime_ext + '+' +  conformation_ext

        ## Add 3' truncation
        if trun3prime != '':
            conformation = conformation + '+' + trun3prime
            conformation_ext = conformation_ext + '+' + trun3prime_ext

        return conformation, conformation_ext

    ## Multiple hits ##
    conformation = ''
    conformation_ext = ''

    for hit_i, hit_j in zip(chain.alignments[:-1], chain.alignments[1:]):

        ## Define orientations
        orientation_i = 'FOR' if hit_i.strand == '+' else 'REV'
        orientation_j = 'FOR' if hit_j.strand == '+' else 'REV'

        ## Compute distance between interval begin positions (i and j) on the consensus == template
        dist_begin = hit_j.tBeg - hit_i.tBeg

        ## Compute distance between interval begin (j) and end (i) on the consensus == template. 
        # Note: these correspond to the bkp junction between both intervals at query level
        dist_junction = hit_j.tBeg - hit_i.tEnd

        #### Define event type connecting both templates
        ## A) Polymerase jump as unexpected order, j before i
        # ---------j----------
        #          ----------i------------
        if dist_begin < 0:
            event = 'JUMP'
            eventLen = dist_begin

        # B) Not jump
        else:

            ## a) Blunt junction
            #  ----------i------------
            #                         ---------j----------
            if (dist_junction == 0):
                event = 'BLUNT'
                eventLen = 'NA'

            ## b) DEL junction
            #  ----------i------------
            #                                ---------j----------
            elif dist_junction >= 1:
                event = 'DEL'
                eventLen = abs(dist_junction)

            ## c) DUP junction
            #  ----------i------------
            #               ---------j----------
            else:
                event = 'DUP'
                eventLen = abs(dist_junction)

        ## Initialize
        if conformation == '':
            conformation = orientation_i + '+' + event + '+' + orientation_j  

            ## Extended conformation
            conformation_ext = orientation_i + '(' + hit_i.alignmentString() + ')+' + event + '(' + str(eventLen) + ')+' + orientation_j + '(' + hit_j.alignmentString() + ')'

        ## Extend       
        else:
            conformation = conformation + '+' + event + '+' + orientation_j

            ## Extended conformation
            conformation_ext = conformation_ext + '+' + event + '(' + str(eventLen) + ')+' + orientation_j + '(' + hit_j.alignmentString() + ')'

    ## Add 5' truncation
    if trun5prime != '':
        conformation = trun5prime + '+' + conformation
        conformation_ext = trun5prime_ext + '+' +  conformation_ext
    
    ## Add 3' truncation
    if trun3prime != '':
        conformation = conformation + '+' + trun3prime
        conformation_ext = conformation_ext + '+' + trun3prime_ext

    return conformation, conformation_ext

        
def compute_perc_resolved(retro_annot):
    '''
    '''
    ## Compute TOTAL polyA length
    if 'POLYA_LEN' in retro_annot:
        polyA_len = sum([int(i) for i in str(retro_annot['POLYA_LEN']).split(',')])
    else:
        polyA_len = 0

    ## a) UNK insertion type
    if (retro_annot['DTYPE_N'] == 'UNK'):
        retro_annot['PERC_RESOLVED'] = 0

    ## b) Solo insertion type
    elif (retro_annot['DTYPE_N'] == 'solo'):

        ## General
        rt_len = int(retro_annot['RT_LEN']) if 'RT_LEN' in retro_annot else 0
        tsd_len = int(retro_annot['TSD_LEN']) if 'TSD_LEN' in retro_annot else 0

        ## Nested
        nested_len = int(retro_annot['NESTED_LEN']) if 'NESTED_LEN' in retro_annot else 0

        ## Hexamer
        hexamer_len = int(retro_annot['HEXAMER_LEN']) if 'HEXAMER_LEN' in retro_annot else 0

        ## Compute templated insertion lengths 
        templated_5prime_len = int(retro_annot['5PRIME_TEMP_LEN']) if '5PRIME_TEMP_LEN' in retro_annot else 0
        templated_3prime_len = int(retro_annot['3PRIME_TEMP_LEN']) if '3PRIME_TEMP_LEN' in retro_annot else 0
        templated_internal_len = int(retro_annot['INTERNAL_TEMP_LEN']) if 'INTERNAL_TEMP_LEN' in retro_annot else 0
        templated_internal_len = int(retro_annot['INTERNAL_TEMP_LEN']) if 'INTERNAL_TEMP_LEN' in retro_annot else 0
        
        ## Compute % resolved
        resolvedLen = tsd_len + rt_len + polyA_len + hexamer_len + nested_len
        resolvedLen = resolvedLen + templated_5prime_len + templated_internal_len + templated_3prime_len
        retro_annot['PERC_RESOLVED'] = resolvedLen / retro_annot['DEL_LEN'] * 100

    ## c) Partnered transduction
    elif (retro_annot['DTYPE_N'] == 'partnered'):

        ## General
        tsd_len = int(retro_annot['TSD_LEN']) if 'TSD_LEN' in retro_annot else 0

        ## Compute RT_LEN:
        if 'CHAIN_repeats' in retro_annot:
            beg, end = retro_annot['CHAIN_repeats'].interval()
            retro_annot['RT_LEN'] = end - beg
        else:
            retro_annot['RT_LEN'] = 0

        ## Hexamer
        hexamer_len = retro_annot['HEXAMER_LEN'] if 'HEXAMER_LEN' in retro_annot else 0

        ## Compute templated insertion lengths 
        templated_5prime_len = int(retro_annot['5PRIME_TEMP_LEN']) if '5PRIME_TEMP_LEN' in retro_annot else 0
        templated_3prime_len = int(retro_annot['3PRIME_TEMP_LEN']) if '3PRIME_TEMP_LEN' in retro_annot else 0
        templated_internal_len = int(retro_annot['INTERNAL_TEMP_LEN']) if 'INTERNAL_TEMP_LEN' in retro_annot else 0

        ## Compute TOTAL 5' transduction length
        if '5PRIME_TD_LEN' in retro_annot:
            td_5prime_len = sum([int(i) for i in retro_annot['5PRIME_TD_LEN'].split(',')])
        else:
            td_5prime_len = 0

        ## Compute TOTAL 3' transduction length
        if '3PRIME_TD_LEN' in retro_annot:
            td_3prime_len = sum([int(i) for i in retro_annot['3PRIME_TD_LEN'].split(',')])
        else:
            td_3prime_len = 0

        ## Compute % resolved
        resolvedLen = td_5prime_len + tsd_len + retro_annot['RT_LEN'] + td_3prime_len + polyA_len + hexamer_len
        resolvedLen = resolvedLen + templated_5prime_len + templated_internal_len + templated_3prime_len

        retro_annot['PERC_RESOLVED'] = resolvedLen / retro_annot['DEL_LEN'] * 100

    ## d) Orphan transduction
    elif (retro_annot['DTYPE_N'] == 'orphan'):
        ## General
        tsd_len = int(retro_annot['TSD_LEN']) if 'TSD_LEN' in retro_annot else 0

        ## Compute % resolved
        resolvedLen = tsd_len + retro_annot['ORPHAN_TD_LEN'] + polyA_len
        retro_annot['PERC_RESOLVED'] = resolvedLen / retro_annot['DEL_LEN'] * 100

    ## e) Processed pseudogene
    elif (retro_annot['DTYPE_N'] == 'PSD'):
        ## General
        tsd_len = int(retro_annot['TSD_LEN']) if 'TSD_LEN' in retro_annot else 0

        ## Compute % resolved
        resolvedLen = tsd_len + retro_annot['PSD_LEN'] + polyA_len
        retro_annot['PERC_RESOLVED'] = resolvedLen / retro_annot['DEL_LEN'] * 100

    ## f) Other insertion type?
    else:
        retro_annot['PERC_RESOLVED'] = 0

    return retro_annot

def solo_annot(retro_annot, chain):
    '''
    '''
    ## Solo if unknown insertion type and family defined
    if (retro_annot['DTYPE_N'] == 'UNK') and ('FAM_N' in retro_annot.keys()):
        retro_annot['DTYPE_N'] = 'solo'

        ## Retrotransposon length
        beg, end = chain.interval()
        retro_annot['RT_LEN'] = end - beg

    return retro_annot

def transductions_5prime(retro_annot, chain_consensus, PAF_ref, exonsBinDb):
    '''
    DO LATER, ONCE fixed issue with SVA alignment into the consensus
    '''
    ## 1. Discard if:
    # - PSD insertions OR
    # - orphan transductions OR
    # - insertions with templated insertions on their 5' ends OR
    # - no hits on the reference OR
    # - Alu insertion OR
    # - not unaligned piece on the consensus longer or equal than Xbp
    minTdLen = 10
    nbTemplates = chain_consensus.nb_templates()[0]
    templates = ','.join(chain_consensus.nb_templates()[1])
    unalignedLen = chain_consensus.interval()[0]

    if (retro_annot['DTYPE_N'] == 'PSD') or (retro_annot['DTYPE_N'] == 'orphan') or ('5PRIME_TEMP_SEQ' in retro_annot) or (len(PAF_ref.alignments) == 0) or ((nbTemplates == 1) and ('Alu' in templates)) or (unalignedLen < minTdLen):
        return retro_annot
    
    ### 2. Search for hits on the reference spanning the non-retrotransposon sequence
    unalignedBeg = 0
    unalignedEnd = unalignedLen
    minPercOverlap = 80
    PAF_overlap = formats.PAF()

    for refHit in PAF_ref.alignments:

        overlap = gRanges.overlap(unalignedBeg, unalignedEnd, refHit.qBeg, refHit.qEnd)[0]

        if not overlap:
            continue

        ## Add hit
        PAF_overlap.alignments.append(refHit)

        ## Check if hit overlaps an exon
        refHit.geneNames = []

        if refHit.tName not in list(exonsBinDb.keys()):
            continue

        exonOverlaps = exonsBinDb[refHit.tName].collect_interval(refHit.tBeg, refHit.tEnd, ['EXON'])

        ## Filter overlaps
        filtered_exonOverlaps = [i for i in exonOverlaps if i[2] >= minPercOverlap]

        ## Collect gene names
        if filtered_exonOverlaps:
            refHit.geneNames = list(set([i[0].optional['geneName'] for i in filtered_exonOverlaps]))

    # 3. Discard as no overlapping hits on the reference
    if not PAF_overlap.alignments:
        return retro_annot

    ## Chain overlapping hits
    chain = PAF_overlap.chain(10, 25)

    ## Collect geneName info
    geneNames = [hit.geneNames for hit in chain.alignments if hit.geneNames]
    geneNames = list(set([i for subList in geneNames for i in subList])) ## Flatten the geneName list
    geneNames = ','.join(geneNames)

    ## Call 5' transduction if enough % assigned
    minPerc = 25
    overlap, overlapLen = gRanges.overlap(unalignedBeg, unalignedEnd, chain.alignments[0].qBeg, chain.alignments[-1].qEnd)
    percAsigned = float(overlapLen) / unalignedLen * 100

    if percAsigned >= minPerc:

        retro_annot['DTYPE_N'] = 'partnered'
        retro_annot['5PRIME_TD_LEN'] = ','.join([str(hit.qEnd - hit.qBeg) for hit in chain.alignments])
        retro_annot['5PRIME_TD_SEQ'] = ','.join([retro_annot['INS_TRIMMED'][hit.qBeg:hit.qEnd] for hit in chain.alignments])
        retro_annot['5PRIME_NB_TD'] = len(chain.alignments)
        retro_annot['5PRIME_TD_COORD'] = ','.join([hit.tName + ':' + str(hit.tBeg) + '-' + str(hit.tEnd) for hit in chain.alignments])
        retro_annot['5PRIME_TD_MAPQ'] = ','.join([str(hit.MAPQ) for hit in chain.alignments])

        beg, end = chain_consensus.interval()
        retro_annot['RT_LEN'] = end - beg

        if geneNames:
            retro_annot['5PRIME_MOB_GENE'] = geneNames

    return retro_annot

def transductions_3prime_hits(retro_annot, chain_consensus, PAF_ref, exonsBinDb):
    '''
    '''
    ## 1. Discard if:
    # - PSD insertions OR
    # - orphan transductions OR
    # - partnered transduction OR
    # - insertions with templated insertions on their 5' ends OR
    # - no hits on the reference OR
    # - Alu insertion OR
    nbTemplates = chain_consensus.nb_templates()[0]
    templates = ','.join(chain_consensus.nb_templates()[1])

    if (retro_annot['DTYPE_N'] == 'PSD') or (retro_annot['DTYPE_N'] == 'orphan') or (retro_annot['DTYPE_N'] == 'partnered') or ('5PRIME_TEMP_SEQ' in retro_annot) or (len(PAF_ref.alignments) == 0) or ((nbTemplates == 1) and ('Alu' in templates)):
        return retro_annot

    ## 2. Discard if:
    # - not unaligned piece on the consensus longer or equal than Xbp
    minTdLen = 10
    seqLen = len(retro_annot['INS_TRIMMED'])
    alignedEnd = chain_consensus.interval()[1]
    unalignedLen = seqLen - alignedEnd

    if (unalignedLen < minTdLen):
        return retro_annot

    ### 3. Search for hits on the reference spanning the non-retrotransposon sequence
    unalignedBeg = alignedEnd
    unalignedEnd = seqLen
    minPercOverlap = 80
    PAF_overlap = formats.PAF()

    for refHit in PAF_ref.alignments:

        overlap = gRanges.overlap(unalignedBeg, unalignedEnd, refHit.qBeg, refHit.qEnd)[0]

        if not overlap:
            continue
        
        ## Add hit
        PAF_overlap.alignments.append(refHit)

        ## Check if hit overlaps an exon
        refHit.geneNames = []

        if refHit.tName not in list(exonsBinDb.keys()):
            continue

        exonOverlaps = exonsBinDb[refHit.tName].collect_interval(refHit.tBeg, refHit.tEnd, ['EXON'])

        ## Filter overlaps
        filtered_exonOverlaps = [i for i in exonOverlaps if i[2] >= minPercOverlap]

        ## Collect gene names
        if filtered_exonOverlaps:
            refHit.geneNames = list(set([i[0].optional['geneName'] for i in filtered_exonOverlaps]))

    # 4. Discard as no overlapping hits on the reference
    if not PAF_overlap.alignments:
        return retro_annot

    ## Chain overlapping hits
    chain = PAF_overlap.chain(10, 25)

    ## Collect geneName info
    geneNames = [hit.geneNames for hit in chain.alignments if hit.geneNames]
    geneNames = list(set([i for subList in geneNames for i in subList])) ## Flatten the geneName list
    geneNames = ','.join(geneNames)

    ## Call 3' transduction if enough % assigned
    minPerc = 25

    overlap, overlapLen = gRanges.overlap(unalignedBeg, unalignedEnd, chain.interval()[0], chain.interval()[1])
    percAsigned = float(overlapLen) / unalignedLen * 100

    if percAsigned >= minPerc:

        retro_annot['DTYPE_N'] = 'partnered'
        retro_annot['3PRIME_TD_LEN'] = ','.join([str(hit.qEnd - hit.qBeg) for hit in chain.alignments])
        retro_annot['3PRIME_TD_SEQ'] = ','.join([retro_annot['INS_TRIMMED'][hit.qBeg:hit.qEnd] for hit in chain.alignments])
        retro_annot['3PRIME_NB_TD'] = len(chain.alignments)
        retro_annot['3PRIME_TD_COORD'] = ','.join([hit.tName + ':' + str(hit.tBeg) + '-' + str(hit.tEnd) for hit in chain.alignments])
        retro_annot['3PRIME_TD_MAPQ'] = ','.join([str(hit.MAPQ) for hit in chain.alignments])

        if geneNames:
            retro_annot['3PRIME_MOB_GENE'] = geneNames

    return retro_annot


def templated_ins(retro_annot, chain):
    '''
    '''
    ## Collect info regarding the templates and hits on the chain
    templates = ','.join(chain.nb_templates()[1])
    nbTemplates = chain.nb_templates()[0]
    nbHits = len(chain.alignments)

    ## Discard inserts aligning in a single template and no hit over retrotransposon consensus database
    if (nbTemplates == 1) and ('consensus' not in templates):
        return retro_annot
    
    ## Define insertion region interval
    offset = 50000
    firstHit = chain.alignments[0]
    insIdRef, coord = firstHit.qName.split(':')
    iRef = insIdRef.split('_')[1]
    iBeg = int(coord) - offset
    iEnd = int(coord) + offset
    
    ## Search for potential templated insertions on the alignment chain
    for idx, hit in enumerate(chain.alignments):
        
        ## Hit adjacent to the insertion region interval
        if (hit.tName == iRef) and (gRanges.overlap(hit.tBeg, hit.tEnd, iBeg, iEnd)[0]):

            ## a) 5' templated
            if idx == 0:
                retro_annot['5PRIME_TEMP_SEQ'] = retro_annot['INS_TRIMMED'][chain.alignments[idx].qBeg:chain.alignments[idx].qEnd]
                retro_annot['5PRIME_TEMP_LEN'] = len(retro_annot['5PRIME_TEMP_SEQ'])
                retro_annot['5PRIME_TEMP_COORD'] = hit.tName + ':' + str(hit.tBeg) + '-' + str(hit.tEnd)
            
            ## b) 3' templated
            elif (idx + 1) == nbHits:
                retro_annot['3PRIME_TEMP_SEQ'] = retro_annot['INS_TRIMMED'][chain.alignments[idx].qBeg:chain.alignments[idx].qEnd]
                retro_annot['3PRIME_TEMP_LEN'] = len(retro_annot['3PRIME_TEMP_SEQ'])
                retro_annot['3PRIME_TEMP_COORD'] = hit.tName + ':' + str(hit.tBeg) + '-' + str(hit.tEnd)

            ## c) Internal templated
            else:
                retro_annot['INTERNAL_TEMP_SEQ'] = retro_annot['INS_TRIMMED'][chain.alignments[idx].qBeg:chain.alignments[idx].qEnd]
                retro_annot['INTERNAL_TEMP_LEN'] = len(retro_annot['INTERNAL_TEMP_SEQ'])
                retro_annot['INTERNAL_TEMP_COORD'] = hit.tName + ':' + str(hit.tBeg) + '-' + str(hit.tEnd)

    return retro_annot

def orphan_td_annot(retro_annot, chain):
    '''
    '''
    templates = ','.join(chain.nb_templates()[1])
    nbTemplates = chain.nb_templates()[0]
    beg, end = chain.interval_template()

    ## Call orphan transduction if:
    # - The insert has no hit on a retrotransposon consensus sequence AND
    # - Insertion type not being already defined
    # - Hits aligning over a single template
    # - Alignment hits encompass at least X% of the insert
    if ('consensus' not in templates) and (retro_annot['DTYPE_N'] == 'UNK') and (nbTemplates == 1) and (chain.perc_query_covered() >= 75):
        retro_annot['DTYPE_N'] = 'orphan'
        retro_annot['ORPHAN_TD_SEQ'] = retro_annot['INS_TRIMMED']
        retro_annot['ORPHAN_TD_LEN'] = len(retro_annot['ORPHAN_TD_SEQ'])
        retro_annot['ORPHAN_TD_MAPQ'] = ','.join([str(hit.MAPQ) for hit in chain.alignments])

        ## TODO: ADD ORPHAN_TD_COORD
        retro_annot['ORPHAN_TD_COORD'] = chain.alignments[0].tName + ':' + str(beg) + '-' + str(end)

        ## Cleanup         
        del retro_annot['INS_TRIMMED']
            
    return retro_annot

def merge_hits(PAFs_A, PAFs_B):
    '''
    '''
    ## Generate a non-redundant list with all insertion ids with available hits
    insIds = list(set(list(PAFs_A.keys()) + list(PAFs_B.keys())))
    PAFs_merged = {}

    for id in insIds:

        ## Initialize PAF object
        PAFs_merged[id] = formats.PAF()

        ## Add A hits if available
        if id in PAFs_A:
            PAFs_merged[id].alignments = PAFs_merged[id].alignments + PAFs_A[id].alignments

        ## Add B hits if available
        if id in PAFs_B:
            PAFs_merged[id].alignments = PAFs_merged[id].alignments + PAFs_B[id].alignments

    return PAFs_merged

def filter_overlapping_hits(inputPAFS, refPAFS, maxPercOverlap):
    '''
    '''
    for insId in inputPAFS.keys():

        if insId not in refPAFS:
            continue

        ## 1. Organize reference hits into a bin database
        refHitsBinDb = hits2binDb(refPAFS[insId].alignments)

        ### 2. Filter input hits redundant to the reference
        inputPAFS[insId].alignments = filter_hits(inputPAFS[insId].alignments, refHitsBinDb, maxPercOverlap)

    return inputPAFS

def hits2binDb(hits):
    '''
    '''
    ## Organize hits into a dictionary 
    hitsDict = {}
    hitsDict['HITS'] = []

    for hit in hits:
        hit.beg = hit.qBeg
        hit.end = hit.qEnd
        hitsDict['HITS'].append(hit)

    ## Create bin database
    binSizes = [100, 1000, 10000]
    binDb = structures.create_bin_database_interval('insert', 0, hits[0].qLen, hitsDict, binSizes)

    return binDb

def filter_hits(hits, refHitsBinDb, maxOverlap):
    '''
    '''
    passHits = []

    for hit in hits:

        overlaps = refHitsBinDb.collect_interval(hit.qBeg, hit.qEnd, ['HITS'])

        if not overlaps:
            passHits.append(hit)
            continue

        overlaps.sort(key=lambda x: x[2], reverse=True)
        percOverlap = overlaps[0][2]

        if percOverlap < maxOverlap:
            passHits.append(hit)
            continue            
    
    return passHits

def alignTrimmedIns2ref(retro_annot, reference, minMAPQ, minHitLen, fileName, outDir):
    '''
    NOTE: Both alignments take ~5'
    '''
    ## 1. Separate long and short insertion events
    retro_annot_short, retro_annot_long = separate_retro_length(retro_annot, 500)

    ## 2. Write fastas with inserted sequences
    FASTA_short_path = outDir + '/' + fileName + '_short.fa'
    fasta_short = trimmed2fasta(retro_annot_short)
    fasta_short.write(FASTA_short_path)

    FASTA_long_path = outDir + '/' + fileName + '_long.fa'
    fasta_long = trimmed2fasta(retro_annot_long)
    fasta_long.write(FASTA_long_path)

    ## 3. Align short inserts with Minimap2   
    preset = 'sr'
    #PAF_short_path = '/Users/brodriguez/Research/Projects/1KGP_LR/SV_annot/resubmission/annot_INS/tmp/trimmed2ref_short.paf'
    PAF_short_path = alignment_minimap2(FASTA_short_path, reference, preset, fileName + '_short', 2, outDir)

    PAF_short_ref = formats.PAF()
    PAF_short_ref.read(PAF_short_path)

    ## 4. Align long inserts with Minimap2   
    preset = 'map-ont'
    #PAF_long_path = '/Users/brodriguez/Research/Projects/1KGP_LR/SV_annot/resubmission/annot_INS/tmp/trimmed2ref_long.paf'
    PAF_long_path = alignment_minimap2(FASTA_long_path, reference, preset, fileName + '_long', 2, outDir)

    PAF_long_ref = formats.PAF()
    PAF_long_ref.read(PAF_long_path)

    ## 5. Merge short and long PAFs
    PAF_ref = formats.PAF()
    PAF_ref.alignments = PAF_short_ref.alignments + PAF_long_ref.alignments

    ## 6. Filter hits based on MAPQ
    PAF_ref.alignments = [hit for hit in PAF_ref.alignments if hit.MAPQ >= minMAPQ]

    ## 7. Filter hits shorter than Xbp
    PAF_ref.alignments = [hit for hit in PAF_ref.alignments if hit.alignmentLen() >= minHitLen]

    ## 8. Generate single PAF objects per inserted sequence:
    PAFs_ref = group_alignments(PAF_ref)

    ## Cleanup
    #unix.rm([FASTA_path, SAM_path, PAF_path, outDir + '/align.err'])

    return PAFs_ref

def alignTrimmedIns2repeats(retro_annot, reference, minMAPQ, minHitLen, fileName, outDir):
    '''
    NOTE: Both alignments take ~5'
    '''
    ## 1. Write fastas with inserted sequences
    FASTA_path = outDir + '/' + fileName + '.fa'
    fasta = trimmed2fasta(retro_annot)
    fasta.write(FASTA_path)

    ## 2. Align trimmed inserts with BWA-mem   
    #PAF_path = '/Users/brodriguez/Research/Projects/1KGP_LR/SV_annot/resubmission/annot_INS/tmp/trimmed2repeats.paf'
    SAM_path = alignment.alignment_bwa(FASTA_path, reference, fileName, 2, outDir)
    PAF_path = alignment.sam2paf(SAM_path, fileName, outDir)

    PAF_repeats = formats.PAF()
    PAF_repeats.read(PAF_path)

    ## 3. Filter hits based on MAPQ
    PAF_repeats.alignments = [hit for hit in PAF_repeats.alignments if hit.MAPQ >= minMAPQ]

    ## 4. Filter hits shorter than Xbp
    PAF_repeats.alignments = [hit for hit in PAF_repeats.alignments if hit.alignmentLen() >= minHitLen]

    ## 5. Generate single PAF objects per inserted sequence:
    PAFs_repeats = group_alignments(PAF_repeats)

    ## Cleanup
    #unix.rm([FASTA_path, SAM_path, PAF_path, outDir + '/align.err'])

    return PAFs_repeats

def alignTrimmedIns2repeats2(retro_annot, reference, minMAPQ, minHitLen, fileName, outDir):
    '''
    NOTE: Both alignments take ~5'
    '''
    ## 1. Write fastas with inserted sequences
    FASTA_path = outDir + '/' + fileName + '.fa'
    fasta = trimmed2fasta(retro_annot)
    fasta.write(FASTA_path)

    ## 2. Align trimmed inserts with BWA-mem  
    #PAF_path = '/Users/brodriguez/Research/Projects/1KGP_LR/SV_annot/resubmission/annot_INS/tmp/trimmed2repeats_lenient.paf' 
    SAM_path = alignment.alignment_bwa(FASTA_path, reference, fileName, 2, outDir)
    PAF_path = alignment.sam2paf(SAM_path, fileName, outDir)

    PAF_repeats = formats.PAF()
    PAF_repeats.read(PAF_path)

    ## 3. Filter hits based on MAPQ
    PAF_repeats.alignments = [hit for hit in PAF_repeats.alignments if hit.MAPQ >= minMAPQ]

    ## 4. Filter hits shorter than Xbp
    PAF_repeats.alignments = [hit for hit in PAF_repeats.alignments if hit.alignmentLen() >= minHitLen]

    ## 5. Generate single PAF objects per inserted sequence:
    PAFs_repeats = group_alignments(PAF_repeats)

    ## Cleanup
    #unix.rm([FASTA_path, SAM_path, PAF_path, outDir + '/align.err'])

    return PAFs_repeats

def nested_annot(retro_annot, reference, minMAPQ, minHitLen, fileName, outDir):
    '''
    NOTE: Both alignments take ~5'
    '''

    ## Create fasta with candidate sequences for nested ins
    FASTA_path = create_fasta_nested_ins(retro_annot, outDir)

    ## Align transduced sequences into the reference
    #PAF_path = '/Users/brodriguez/Research/Projects/1KGP_LR/SV_annot/resubmission/annot_INS/tmp/td3prime2consensus.paf'
    SAM_path = alignment.alignment_bwa(FASTA_path, reference, 'td3prime2consensus', 2, outDir)     
    PAF_path = alignment.sam2paf(SAM_path, 'td3prime2consensus', outDir)

    ## Read alignments and create chain
    PAF_cons = formats.PAF()
    PAF_cons.read(PAF_path)

    ## Filter out hits with MAPQ == 0
    PAF_cons.alignments = [hit for hit in PAF_cons.alignments if hit.MAPQ > 0]
    
    ## Discard if no hit passes the filters
    if not PAF_cons.alignments:
        return
    
    ## Group alignments by insertion id
    consHitsTd = group_alignments(PAF_cons)

    ## For each insertion with nested sequence aligning on a consensus
    for insId in consHitsTd:
 
        if 'FAM_N' not in retro_annot[insId]:
            continue

        ## create chain
        chain = consHitsTd[insId].chain(10, 25)

        ## annotate conformation
        conf, conf_ext = annotate_conformation(chain)

        ## select last operation in the conformation string
        last_op = conf.split('+')[-1]
        families = list(set([hit.tName.split('|')[1] for hit in chain.alignments]))
        family = ','.join(families)

        if (chain.perc_query_covered() > 80) and (last_op != 'TRUN'):

            ## Create nested entry
            retro_annot[insId]['DTYPE_N'] = 'solo'
            retro_annot[insId]['NESTED_LEN'] = retro_annot[insId]['3PRIME_TD_LEN']
            retro_annot[insId]['NESTED_SEQ'] = retro_annot[insId]['3PRIME_TD_SEQ']
            retro_annot[insId]['NESTED_FAM'] = family
            retro_annot[insId]['NESTED_COORD'] = ','.join([hit.tName + ':' + str(hit.tBeg) + '-' + str(hit.tEnd) for hit in chain.alignments])
            retro_annot[insId]['NESTED_MAPQ'] = ','.join([str(hit.MAPQ) for hit in chain.alignments]) 

            if 'CHAIN_repeats' in retro_annot[insId]:
                beg, end = retro_annot[insId]['CHAIN_repeats'].interval()
                retro_annot[insId]['RT_LEN'] = end - beg
            else:
                retro_annot[insId]['RT_LEN'] = 0

            ## Remove 3' transduction info
            retro_annot[insId].pop('3PRIME_TD_LEN', None)
            retro_annot[insId].pop('3PRIME_TD_SEQ', None)
            retro_annot[insId].pop('3PRIME_NB_TD', None)
            retro_annot[insId].pop('3PRIME_TD_COORD', None)
            retro_annot[insId].pop('3PRIME_TD_MAPQ', None)

    return retro_annot

def separate_retro_length(retro_annot, minLen):
    '''
    '''
    ## Initialize dict
    retro_annot_short = {}
    retro_annot_long = {}

    ## Sparate events by length
    for id, annot in retro_annot.items():
        length = len(annot['INS_TRIMMED'])

        ## Short 
        if length <= minLen:
            retro_annot_short[id] = annot
        ## Long
        else:
            retro_annot_long[id] = annot
    
    return retro_annot_short, retro_annot_long


def family_annot_insert(retro_annot, chain):
    '''
    '''
    ## 1. Discard processed pseudogenes and orphan transductions
    if (retro_annot['DTYPE_N'] == 'PSD') or (retro_annot['DTYPE_N'] == 'orphan'):
        return retro_annot

    ## 2. Determine insertion family
    families = list(set([hit.tName.split('|')[1] for hit in chain.alignments]))
    nbFam = len(families)

    ### a) SVA 
    if ('SVA' in families):
        retro_annot['FAM_N'] = 'SVA'

    ### b) Non SVA retroelement
    elif nbFam == 1:
        retro_annot['FAM_N'] = families[0]

    ### c) Chimeric retroelement insertion
    elif (nbFam > 1):
        retro_annot['DTYPE_N'] = 'chimera'
        retro_annot['FAM_N'] = ','.join(families)

    ### d) Unresolved family

    return retro_annot

def ins2fasta(vcf):
    '''
    Create fasta object containing insertions in a VCF file

    Input:
        1. vcf: Path to VCF file

    Output: 
        1. fasta object
    '''     
    ## 1. Initialize fasta object
    fasta = formats.FASTA() 

    ## 2. Collect inserted sequences
    for variant in vcf.variants:

        ## Exclude insertions longer than 20Kb
        if len(variant.ref) > 20000:
            continue
    
        ## Save into dictionary
        insId = variant.ID + '_' + variant.chrom + ':' + str(variant.pos) 
        seq = variant.ref
        fasta.seqDict[insId] = seq

    return fasta

def trimmed2fasta(annot):
    '''
    Organize inserted sequences into a fasta object

    Input:
        1. annot: 
        2. outDir: Output directory

    Output: 
        1. fasta: path to fasta object inserted sequences
    '''     
    ## 1. Initialize fasta object
    fasta = formats.FASTA() 

    ## 2. Collect inserted sequences
    for insId in annot.keys():
        fasta.seqDict[insId] = annot[insId]['INS_TRIMMED']

    return fasta

def group_alignments(paf):
    '''
    '''     
    pafDict = {}

    ## For each hit
    for hit in paf.alignments:

        # Initialize paf object for this inserted sequence
        if hit.qName not in pafDict:
            pafDict[hit.qName] = formats.PAF()
    
        # Add hit to the corresponding paf
        pafDict[hit.qName].alignments.append(hit)

    return pafDict

def align_ins2ref(VCF, reference, outDir):
    '''
    '''
    ## Separate long and short insertion events
    short_VCF, long_VCF = separate_events_length(VCF, 500)

    ## Write fastas with inserted sequences
    FASTA_short_path = outDir + '/inserted_sequences_short.fa'
    fasta_short = ins2fasta(short_VCF)
    fasta_short.write(FASTA_short_path)

    FASTA_long_path = outDir + '/inserted_sequences_long.fa'
    fasta_long = ins2fasta(long_VCF)
    fasta_long.write(FASTA_long_path)

    ## Align short inserts with Minimap2   
    #PAF_short_path = '/Users/brodriguez/Research/Projects/1KGP_LR/SV_annot/resubmission/annot_INS/tmp/inserts2ref_short.paf'
    preset = 'sr'
    PAF_short_path = alignment_minimap2(FASTA_short_path, reference, preset, 'inserts2ref_short', 2, outDir)

    PAF_short_ref = formats.PAF()
    PAF_short_ref.read(PAF_short_path)

    ## Align long inserts with Minimap2   
    #PAF_long_path = '/Users/brodriguez/Research/Projects/1KGP_LR/SV_annot/resubmission/annot_INS/tmp/inserts2ref_long.paf'
    preset = 'map-ont'
    PAF_long_path = alignment_minimap2(FASTA_long_path, reference, preset, 'inserts2ref_long', 2, outDir)

    PAF_long_ref = formats.PAF()
    PAF_long_ref.read(PAF_long_path)

    ## Merge short and long PAFs
    PAF_ref = formats.PAF()
    PAF_ref.alignments = PAF_short_ref.alignments + PAF_long_ref.alignments

    ## Group alignments by insertion
    refHitsIns = group_alignments(PAF_ref)

    return refHitsIns
    
def alignment_minimap2(INPUT, reference, preset, fileName, processes, outDir):
    '''
    '''
    PAF = outDir + '/' + fileName + '.paf'
    err = open(outDir + '/' + fileName + '.err', 'w') 
    command = 'minimap2 -t ' + str(processes) + ' -x ' + preset + ' ' + reference + ' ' + INPUT + ' > ' + PAF
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Local alignment failed' 
        log.step(step, msg)

    return PAF

def separate_events_length(VCF, maxLen):
    '''
    '''

    ### Create VCFs wiht short and long events
    short_VCF = formats.VCF()
    short_VCF.header = VCF.header

    long_VCF = formats.VCF()
    long_VCF.header = VCF.header

    ### Separate events by size
    for variant in VCF.variants:

        ## Compute variant len
        length = len(variant.ref)

        ## Short events
        if length <= maxLen:
            short_VCF.variants.append(variant)

        ## Length events
        else:
            long_VCF.variants.append(variant)

    return short_VCF, long_VCF
    
def filter_VCF_len(VCF, maxLen):
    '''
    '''
    ### Create filtered VCF
    filtered_VCF = formats.VCF()
    filtered_VCF.header = VCF.header
 
    ### Apply filtering
    for variant in VCF.variants:

        ## Compute variant len
        length = len(variant.ref)

        ## Filter out variant if len > threshold 
        if length <= maxLen:
            filtered_VCF.variants.append(variant)
    
    return filtered_VCF

######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('vcf', help='Path to VCF file with assembly-based SV calls')
parser.add_argument('TRF', help='Path to Tandem Repeat Finder (TRF) output file for the insertions')
parser.add_argument('vntrDb', help='Path to database of annotated VNTR loci')
parser.add_argument('exons', help='Path to bed file containing exon annotations')
parser.add_argument('repeatsAnnot', help='Path to bed file containing repetitive element annotations')
parser.add_argument('repeatsDb', help='Path to FASTA file containing sequences for repeats')
parser.add_argument('reference', help='Path to FASTA file containing the reference genome')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('--replace-info', action="store_true", default=False, dest='replaceInfo', help='Replace info fields by RETRO-related ones. Default: False (== add RETRO to pre-existing fields)')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory' )

## 2. Parse user input ##
args = parser.parse_args()
vcf = args.vcf
TRF_out = args.TRF
vntrDb = args.vntrDb
exons = args.exons
repeatsAnnot = args.repeatsAnnot
repeatsDb = args.repeatsDb
reference = args.reference
fileName = args.fileName
replaceInfo = args.replaceInfo
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('vcf: ', vcf)
print('TRF: ', TRF_out)
print('vntrDb: ', vntrDb)
print('exons: ', exons)
print('repeatsAnnot: ', repeatsAnnot)
print('repeatsDb: ', repeatsDb)
print('reference: ', reference)
print('fileName: ', fileName)
print('replaceInfo: ', replaceInfo)
print('outDir: ', outDir, "\n")


##########
## MAIN ##
##########

## 0.A Create temporary directory
#################################
print('0.A Create temporary directory')

tmpDir = outDir + '/tmp'
unix.mkdir(tmpDir)

## 0.B Add deletion length to each variant
###########################################
### Read VCF (6" with no genotypes)
VCF = formats.VCF()
VCF.read(vcf)

### Update variant info field by adding DEL_LEN 
for variant in VCF.variants:
    annot = {}
    annot['DEL_LEN'] = len(variant.ref) 
    variant.info.update(annot)

## 0.C Filter out DEL longer than a threshold
##############################################
### VCF filtering
maxLen = 50000

filtered_VCF = filter_VCF_len(VCF, maxLen)
print('VCF: ', len(VCF.variants), len(filtered_VCF.variants))

## 1. Load input files 
#######################
print('1. Load input files')

### Create bin database containing annotated VNTRs (14")
print('1.a VNTR database')
#vntrBinDb = formats.bed2binDb(vntrDb, filtered_VCF.header.refLengths, 2)
vntrBinDb = None

### Create bin database containing exons (30")
print('1.b Annotated exons')
exonsBinDb = formats.bed2binDb(exons, filtered_VCF.header.refLengths, 2)

## 2. Align inserts to reference
#################################
print('2. Align inserts to reference')
PAFs_ref = align_ins2ref(filtered_VCF, reference, tmpDir)

## 3. Annotate retrotransposed insertions
##########################################
print('3. Annotate retrotransposed insertions')
retro_VCF, unk_VCF = annot_retro(filtered_VCF, exonsBinDb, PAFs_ref, reference, repeatsDb, TRF_out, tmpDir)

print('NB_retro_VCF: ', len(retro_VCF.variants))
print('NB_UNK: ', len(unk_VCF.variants))

## 4. Annotate VNTR events
###########################
print('4. Annotate VNTR events')
vntr_VCF, unk_VCF = annot_vntr(unk_VCF, TRF_out, vntrBinDb)
print('NB_VNTR: ', len(vntr_VCF.variants))
print('NB_UNK: ', len(unk_VCF.variants))

## 5. Print UNK
##################

## For each variant assess if numt
#for variant in unk_VCF.variants:
 
    #insId = variant.ID + '_' + variant.chrom + ':' + str(variant.pos) 
    #print('REMAINING_UNK: ', insId, len(variant.ref), variant.ref)

## 6. Write output VCF files with annotations
##############################################
print('6. Write output VCF files with annotations')

## Create output VCF
out_VCF = formats.VCF()

## Add header
out_VCF.header = VCF.header

## Add annotated variants
out_VCF.variants = unk_VCF.variants + retro_VCF.variants + vntr_VCF.variants

## Define RETRO annotation info fields
info2add = {'DTYPE_N': ['.', 'String', 'Type of insertion: solo, partnered, orphan, PSD (Processed pseudogene), VNTR'], \
            'DEL_LEN': ['1', 'Integer', 'Insertion length'], \
            'STRAND': ['.', 'String', 'Insertion DNA strand (+ or -)'], \
            'PERC_RESOLVED': ['1', 'Float', 'Percentage of the inserted sequence that has been resolved'], \
            'FAM_N': ['.', 'String', 'Transposon family'], \
            'RT_LEN': ['1', 'Integer', 'Length of the piece of the inserted sequence corresponding to a transposon'], \
            
            ## Annotated repeat at the insertion breakpoint
            'REPEAT_BKP': ['.', 'String', 'Reference repeat at the insertion breakpoint. Format: class,family,repBegCoord,repEndCoord,strand,milliDiv'], \

            ## Conformation and structure related fields
            'CONFORMATION': ['.', 'String', 'Insertion conformation'], \
            'CONFORMATION_EXT': ['.', 'String', 'Detailed insertion conformation, including the length for each component'], \
            'NOT_CANONICAL': ['0', 'Flag', 'Not canonical MEI'], \

            ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">

            ## TSD specific fields
            'TSD_LEN': ['1', 'Integer', 'Target site duplication length'], \
            'TSD_SEQ': ['.', 'String', 'Target site duplication sequence'], \
            
            ## PolyA specific fields
            'POLYA_LEN': ['1', 'Integer', 'Poly(A) length'], \
            'POLYA_SEQ': ['.', 'String', 'Poly(A) sequence'], \

            ## PSD specific fields
            'SRC_GENE': ['.', 'String', 'Source Gene Name'], \
            'NB_EXONS': ['1', 'Integer', 'Number of exons retrotransposed, processed pseudogene'], \

            ## Source elements specific fields
            'SOURCE_COORD': ['.', 'String', 'Source element genomic coordinates chr:beg-end'], \
            'SOURCE_CYTOID': ['.', 'String', 'Source element cytoband identifier'], \
            
            ## Orphan transductions specific fields
            'ORPHAN_TD_LEN': ['1', 'Integer', 'Orphan transduction length'], \
            'ORPHAN_TD_COORD': ['.', 'String', 'Orphan transduction coordinates'], \
            'ORPHAN_TD_MAPQ': ['.', 'String', 'MAPQ for orphan transduction hit'], \
            'ORPHAN_TD_SEQ': ['.', 'String', 'Orphan transduction sequence'], \
            
            ## 5' transductions
            '5PRIME_NB_TD': ['1', 'Integer', 'Number of 5-prime transductions'], \
            '5PRIME_TD_LEN': ['1', 'Integer', 'Lengths for 5-prime transductions'], \
            '5PRIME_TD_COORD': ['.', 'String', 'Genomic coordinates for 5-prime transductions'], \
            '5PRIME_TD_MAPQ': ['.', 'String', 'MAPQ for 5-prime transduction hits'], \
            '5PRIME_TD_SEQ': ['.', 'String', 'Sequences for 5-prime transductions'], \
            '5PRIME_MOB_GENE': ['.', 'String', 'Gene name that is mobilized via a 5-prime transduction'], \

            ## 3' transductions
            '3PRIME_NB_TD': ['1', 'Integer', 'Number of 3-prime transductions'], \
            '3PRIME_TD_LEN': ['1', 'Integer', 'Lengths for 3-prime transductions'], \
            '3PRIME_TD_COORD': ['.', 'String', 'Genomic coordinates for 3-prime transductions'], \
            '3PRIME_TD_MAPQ': ['.', 'String', 'MAPQ for 3-prime transduction hits'], \
            '3PRIME_TD_SEQ': ['.', 'String', 'Sequences for 3-prime transductions'], \
            '3PRIME_MOB_GENE': ['.', 'String', 'Gene name that is mobilized via a 3-prime transduction'], \
            
            ## 5' templated insertion
            '5PRIME_TEMP_LEN': ['1', 'Integer', 'Length for 5-prime templated insertion'], \
            '5PRIME_TEMP_COORD': ['.', 'String', 'Genomic coordinates for 5-prime templated insertion'], \
            '5PRIME_TEMP_SEQ': ['.', 'String', 'Sequence for 5-prime templated insertion'], \

            ## Internal templated insertion
            'INTERNAL_TEMP_LEN': ['1', 'Integer', 'Length for internal templated insertion'], \
            'INTERNAL_TEMP_COORD': ['.', 'String', 'Genomic coordinates for internal templated insertion'], \
            'INTERNAL_TEMP_SEQ': ['.', 'String', 'Sequence for internal templated insertion'], \

            ## 3' templated insertion
            '3PRIME_TEMP_LEN': ['1', 'Integer', 'Length for 3-prime templated insertion'], \
            '3PRIME_TEMP_COORD': ['.', 'String', 'Genomic coordinates for 3-prime templated insertion'], \
            '3PRIME_TEMP_SEQ': ['.', 'String', 'Sequence for 3-prime templated insertion'], \

            ## Nested insertion event
            'NESTED_LEN': ['1', 'Integer', 'Lengths for nested 3-prime mei event'], \
            'NESTED_FAM': ['.', 'String', 'Family for nested mei'], \
            'NESTED_COORD': ['.', 'String', 'MEI consensus coordinates for nested mei'], \
            'NESTED_MAPQ': ['.', 'String', 'MAPQ for nested mei hits'], \
            'NESTED_SEQ': ['.', 'String', 'Sequence for nested mei'], \

            ## SVA specific fields
            'HEXAMER_SEQ': ['.', 'String', 'Hexamer sequence at the 5-prime end'], \
            'HEXAMER_LEN': ['1', 'Integer', 'Hexamer sequence length'], \
            'VNTR_LEN': ['1', 'Integer', 'VNTR sequence length'], \
            'VNTR_COORD': ['.', 'String', 'VNTR begin and end position in the insert'], \
            
            ## VNTR specific fields
            'MOTIFS': ['.', 'String', 'Comma separated list of motifs within a VNTR'], \
            'NB_MOTIFS': ['1', 'Integer', 'Number of motifs within a VNTR'], \
            
            ## NUMT specific fields
            'MT_COORD': ['.', 'String', 'Mitochrondrial coordinates the inserted sequence derives from'], \
            
            ## DUP specific fields
            'DUP_COORD': ['.', 'String', 'Genomic coordinates for the duplicated sequence']
            }

info2add_order1 = []
info2add_order2 = ['DTYPE_N', 'DEL_LEN', 'PERC_RESOLVED', 'MT_COORD', 'NB_MOTIFS', 'MOTIFS', 'DUP_COORD', 'STRAND', 'FAM_N', 'REPEAT_BKP', 'NOT_CANONICAL', 'CONFORMATION', 'CONFORMATION_EXT', 'RT_LEN', 'TSD_LEN', 'TSD_SEQ', 'HEXAMER_LEN', 'HEXAMER_SEQ', 'VNTR_LEN', 'VNTR_COORD', 'POLYA_LEN', 'POLYA_SEQ', 'SRC_GENE', 'NB_EXONS', 'SOURCE_COORD', 'SOURCE_CYTOID', 'ORPHAN_TD_LEN', 'ORPHAN_TD_COORD', 'ORPHAN_TD_MAPQ', 'ORPHAN_TD_SEQ', '5PRIME_NB_TD', '5PRIME_TD_LEN', '5PRIME_TD_COORD', '5PRIME_TD_MAPQ', '5PRIME_TD_SEQ', '5PRIME_MOB_GENE', '3PRIME_NB_TD', '3PRIME_TD_LEN', '3PRIME_TD_COORD', '3PRIME_TD_MAPQ', '3PRIME_TD_SEQ', '3PRIME_MOB_GENE', '5PRIME_TEMP_LEN', '5PRIME_TEMP_COORD', '5PRIME_TEMP_SEQ', 'INTERNAL_TEMP_LEN', 'INTERNAL_TEMP_COORD', 'INTERNAL_TEMP_SEQ', '3PRIME_TEMP_LEN', '3PRIME_TEMP_COORD', '3PRIME_TEMP_SEQ', 'NESTED_LEN', 'NESTED_FAM', 'NESTED_COORD', 'NESTED_MAPQ', 'NESTED_SEQ']

for i in out_VCF.info_order:
    if i not in info2add_order2:
        info2add_order1.append(i)

## Add info fields with RETRO annotations to the VCF header
if replaceInfo:
    ## Replace info
    out_VCF.header.info = info2add

    ## Replace order list
    out_VCF.info_order = info2add_order2

else:
    ## Append info
    out_VCF.header.info.update(info2add)

    ## Append to the order list
    out_VCF.info_order = info2add_order1 + info2add_order2

## Write VCF file
out_VCF.write(out_VCF.info_order, out_VCF.format_order, fileName , outDir)

## Cleanup
#unix.rm([tmpDir])
