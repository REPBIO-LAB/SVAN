## DEPENDENCIES ##
# External
import os
import pandas as pd 
import math 

# Internal
from GAPI import unix
from GAPI import formats
from GAPI import alignment
from GAPI import bamtools

def snp2haplotype(vcf, lenght):
    '''
    Generate haplotype vector from a VCF containing SNPs calls.

    Input: 
        1. vcf: SNP and INDEL calls in VCF format
        2. length: length of reference sequence

    Output:
        1. haplotype: haplotype vector containing the status for each reference position 
                      (0, no change; 1-12: one of the 12 possible single base substitutions)
    '''
    ## 1. Create dictionary with substitution type codes:
    codes = {
            'A->T': 1,
            'A->G': 2,
            'A->C': 3,
            'T->A': 4,
            'T->G': 5,
            'T->C': 6,
            'G->A': 7,
            'G->T': 8,
            'G->C': 9,
            'C->A': 10,
            'C->T': 11,
            'C->G': 12           
            }

    ## 2. Initialize haplotype vector as always reference (0s)
    haplotype = [0] * lenght

    ## 3. Update haplotype vector with identified SNPs
    # Read VCF
    VCF = formats.VCF()
    VCF.read(vcf)

    # Iterate over each identified variant
    for variant in VCF.variants:

        ## Skip reference calls and indels
        if variant.alt == '.' or 'INDEL' in variant.info:
            continue
    
        ## Add SNP call to the vector
        key = variant.ref + '->' + variant.alt
    
        ## Skip if unknown nucleotyde substitution or multiallelic
        if key not in codes:
            continue

        code = codes[key]
        index = variant.pos - 1 # Convert 1-based (vcf) to 0-based (vector)
        haplotype[index] = code

    return haplotype

def indel2haplotype(vcf, lenght):
    '''
    Generate haplotype vector from a VCF containing INDEL calls.

    Input: 
        1. vcf: SNP and INDEL calls in VCF format
        2. length: length of reference sequence

    Output:
        1. haplotype: haplotype vector containing the status for each reference position 
                      (0, no change; -X: deletion; +X: insertion)
    '''
    ## 2. Initialize haplotype vector as always reference (0s)
    haplotype = [0] * lenght

    ## 3. Update haplotype vector with identified SNPs
    # Read VCF
    VCF = formats.VCF()
    VCF.read(vcf)

    # Iterate over each identified variant
    for variant in VCF.variants:

        # a) Skip if not indel
        if 'INDEL' not in variant.info:
            continue

        # a) Deletion
        elif len(variant.ref) > len(variant.alt):

            index = variant.pos - 1 # Convert 1-based (vcf) to 0-based (vector)
            haplotype[index] = '-' + variant.ref

        # b) Insertion
        elif len(variant.ref) < len(variant.alt):
            index = variant.pos - 1 # Convert 1-based (vcf) to 0-based (vector)
            haplotype[index] = '+' + variant.alt

    return haplotype


def compute_haplotypes(sequences, reference, outDir):
    '''
    Compute haplotype matrix for a set of input sequences by comparing against a reference sequence
    
    Input:
        1. sequences: path to FASTA file containing input sequences
        2. reference: path to FASTA file containing reference sequence
        3. outDir: output file directory
    '''
    ## 1. Create output firectory ##
    unix.mkdir(outDir)

    ## 2. Compute reference length ##
    refFasta = formats.FASTA()             
    refFasta.read(reference)
    refLen = len(list(refFasta.seqDict.values())[0])

    ## 3. Compute haplotypes ##
    ## Read input sequences
    inputFasta = formats.FASTA()             
    inputFasta.read(sequences)

    snp_haplotypes = {}
    indel_haplotypes = {}

    ## For each input sequence
    for seqId, seq in inputFasta.seqDict.items():
        
        ## 1. Create tmp output directory
        tmpDir = outDir + '/' + seqId
        unix.mkdir(tmpDir)

        ## 2. Create fasta containing target sequence
        targetFasta = formats.FASTA()             
        targetFasta.seqDict[seqId] = seq
        target = tmpDir + '/' + seqId + '.fa'
        targetFasta.write(target)

        ## 3. Map target sequence against consensus
        SAM = alignment.alignment_bwa(target, reference, seqId, 1, tmpDir)

        ## 4. SAM to BAM conversion and sorting
        BAM = bamtools.SAM2BAM(SAM, tmpDir)

        ## 5. Call SNPs and INDELS
        calls = tmpDir + '/' + seqId + '.vcf' 
        command = 'samtools mpileup -f ' + reference + ' -g ' + BAM + ' | bcftools call --ploidy 1 -c - > ' + calls
        os.system(command) 
        
        ## 6. Convert SNP calls into vector representation
        snp_haplotype = snp2haplotype(calls, refLen)
        snp_haplotypes[seqId] = snp_haplotype

        ## 7. Convert INDEL calls into vector representation
        indel_haplotype = indel2haplotype(calls, refLen)
        indel_haplotypes[seqId] = indel_haplotype

        ## Cleanup
        unix.rm([tmpDir])

    ## Generate dataframe containing haplotypes
    colNames = list(str(i) for i in range(1, refLen + 1))
    snp_haplotypesDf = pd.DataFrame.from_dict(snp_haplotypes, orient='index', columns=colNames)
    indel_haplotypesDf = pd.DataFrame.from_dict(indel_haplotypes, orient='index', columns=colNames)

    ## Cleanup
    #unix.rm([outDir])

    return snp_haplotypesDf, indel_haplotypesDf


def compute_haplotypes_minimap(sequences, ref, index, outDir):
    '''
    Compute haplotype matrix for a set of input sequences by comparing against a reference sequence
    
    Input:
        1. sequences: path to FASTA file containing input sequences
        2. ref: path to FASTA file containing reference
        3. index: path to minimap2 index for reference 
        4. outDir: output file directory
    '''
    ## 1. Create output firectory ##
    unix.mkdir(outDir)

    ## 2. Compute reference length ##
    fasta_ref = formats.FASTA()             
    fasta_ref.read(ref)
    refLen = len(list(fasta_ref.seqDict.values())[0])

    ## 3. Compute haplotypes ##
    ## Read input sequences
    inputFasta = formats.FASTA()             
    inputFasta.read(sequences)

    snp_haplotypes = {}
    indel_haplotypes = {}

    ## For each input sequence
    for seqId, seq in inputFasta.seqDict.items():
        
        ## 3.1 Create tmp output directory
        tmpDir = outDir + '/' + seqId
        unix.mkdir(tmpDir)

        ## 3.2 Create fasta containing target sequence
        targetFasta = formats.FASTA()             
        targetFasta.seqDict[seqId] = seq
        target = tmpDir + '/' + seqId + '.fa'
        targetFasta.write(target)

        ## 3.3 Map target sequence against consensus
        BAM = alignment.alignment_minimap2_bam(target, index, seqId, 1, tmpDir)

        ## 3.4 Call SNPs and INDELS
        calls = tmpDir + '/' + seqId + '.vcf' 
        command = 'samtools mpileup -f ' + ref + ' -g ' + BAM + ' | bcftools call --ploidy 1 -c - > ' + calls
        os.system(command) 
        
        ## 3.5 Convert SNP calls into vector representation
        snp_haplotype = snp2haplotype(calls, refLen)
        snp_haplotypes[seqId] = snp_haplotype

        ## 3.6 Convert INDEL calls into vector representation
        indel_haplotype = indel2haplotype(calls, refLen)
        indel_haplotypes[seqId] = indel_haplotype

        ## Cleanup
        unix.rm([tmpDir])

    ## 4. Generate dataframe containing haplotypes
    colNames = list(str(i) for i in range(1, refLen + 1))
    snp_haplotypesDf = pd.DataFrame.from_dict(snp_haplotypes, orient='index', columns=colNames)
    indel_haplotypesDf = pd.DataFrame.from_dict(indel_haplotypes, orient='index', columns=colNames)


    ## 5. Cleanup
    unix.rm([outDir])

    return snp_haplotypesDf, indel_haplotypesDf

def diagnostic_scores(haplotypes):
    '''
    '''
    scores = pd.DataFrame(index=haplotypes.index, columns=haplotypes.columns)

    ## Iterate over each position on the haplotypes
    for position in haplotypes:
    
        counts_norm = haplotypes[position].value_counts(normalize=True)
        score_map = counts_norm.apply(lambda x: 1 - x)

        scores_pos = pd.Series([]) 

        for index, value in haplotypes[position].items():
            scores_pos[index] = score_map[value]

        scores[position] = scores_pos

    return scores

def compare_haplotypes(targetHaplo, refHaplo, haploScores):
    '''
    '''
    maxScore = haploScores.sum()
    score = 0

    for pos, status in targetHaplo.items():

        if pos in refHaplo:
            refStatus = refHaplo[pos]
            posScore = haploScores[pos]

            if status == refStatus:
                score += posScore
            
            else:
                score -= posScore
    
    ## Normalize score by the maximum possible
    normScore = score / maxScore

    ## Map score from -1 - 1 space to 0 - 1 space
    normScore = (normScore + 1) / 2
    normScore = round(normScore, 2)

    return normScore


def assign_haplotypes(targetHaplos, refHaplos, refScores):
    '''
    '''
    assignations = pd.DataFrame(index=targetHaplos.index, columns=refHaplos.index)
    cols = list(refHaplos.index)

    ## For each target haplotype
    for targetId, targetHaplo in targetHaplos.iterrows():

        ## For each reference haplotype:
        for refId, refHaplo in refHaplos.iterrows():

            ## Select the score array corresponding to the reference haplotype
            haploScores = refScores.loc[refId, :]

            ## Compare target and reference haplotypes
            assignations.loc[targetId, refId] = compare_haplotypes(targetHaplo, refHaplo, haploScores)

    ## Select max scoring source (if several possible set as None)
    for targetId, scores in assignations.iterrows():

        ## Initialize as None
        assignations.loc[targetId, 'first'] = 'None'
        assignations.loc[targetId, 'lh1'] = 'None'
        assignations.loc[targetId, 'second'] = 'None'
        assignations.loc[targetId, 'lh2'] = 'None'            
        assignations.loc[targetId, 'logRatio'] = 'None'  

        ## Select maximum scoring source elements
        maxScore = scores.max()
        maxSources = scores[scores == maxScore]

        ## Skip assignation if multiple maximum scoring source elements
        if maxSources.size == 1:
        
            ## Select second scoring source
            secondScore = scores[scores != maxScore].max()

            ## Compute the log2 ratio between max and second likely source
            ratio = maxScore / secondScore 
            logRatio = round(math.log(ratio, 2), 2)

            ## Skip assignation if log ratio < 0.08
            if logRatio >= 0.08:
                assignations.loc[targetId, 'first'] = list(maxSources.index)[0]
                assignations.loc[targetId, 'lh1'] = maxScore

                assignations.loc[targetId, 'second'] = list(scores[scores == secondScore].index)[0]
                assignations.loc[targetId, 'lh2'] = secondScore

                assignations.loc[targetId, 'logRatio'] = logRatio

    ## Reorder column names
    ordered = ['first', 'lh1', 'second', 'lh2', 'logRatio'] + cols
    assignations = assignations[ordered] 

    return assignations
    
    
def remove_noninformative(haplo, scores):
    '''
    '''
    ## Drop non-informative positions at score matrix
    filtered_scores = scores.loc[:, (scores != 0).any(axis=0)]

    ## Filter haplotype matrix accordingly
    filtered_haplo = haplo[filtered_scores.columns]

    return filtered_haplo, filtered_scores