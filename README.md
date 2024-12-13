## Structural Variants ANnotator (SVAN)

SVAN is a computational method for the annotation and classification of sequence-resolved insertions and deletions based on their sequence features into distinct classes, including Mobile Element Insertions (MEI), processed pseudogene integrations, various forms of duplications, tandem repeats expansions/contractions and nuclear-mitochondrial segments (NUMT). It primarily takes a VCF file containing the insertion or deletion SV calls as input and produces a second VCF with the annotations for each SV. It is compatible with the output of any long-read SV caller, as long as the sequence for the insertion or deletion events are included in the "alt" and "ref" field of the input VCF, respectively. 
 
SVAN has been used for the SV characterization of 1,019 samples sequenced with long-reads from the 1000 Genomes Project:

Schloissnig et al., “Long-read sequencing and structural variant characterization in 1,019 samples from the 1000 Genomes Project”, bioRxiv, April 20, 2024, https://doi.org/10.1101/2024.04.18.590093

## Download 
Two different ways:

* Go to the releases tab and download the latest release. 

* Clone the git repository in case you want the latest version of the code:

```
# Move into the folder in which you want to clone the repositoy.
$ cd ~/apps
# Clone it.
$ git clone https://github.com/REPBIO-LAB/SVAN.git 
```

SVAN does not require any further installation step. It is written in Python and can be run as a standalone application on diverse Linux systems. 

## Requirements
1. Hardware:

    * 64-bits CPU

2. Software:

    * 64-bit Linux System
    * Python v3.5.4 or higher
    * paftools v.r755 (https://github.com/lh3/minimap2/tree/master/misc)
    * bwa-mem v0.7.17-r1188 (https://github.com/lh3/bwa)
    * minimap2 v2.10-r764-dirty (https://github.com/lh3/minimap2)

3. Python libraries 
    * pysam 
    * cigar
    * itertools
    * Biopython
    * subprocess
    * pandas
    * scipy
    * numpy

## Input
SVAN takes as input 6 mandatory arguments:

   1. VCF: Input VCF file containing sequence-resolved insertion or deletion SV calls. 
   2. TRF: Output for Tandem Repeat Finder (TRF) execution on the inserted or deleted sequence for each SV in the input VCF 
   3. VNTR: Bed file containing VNTR annotation on the reference
   4. EXONS: Bed file containing EXON annotation on the reference
   5. REPEATS: Bed file containing REPEATS annotated with RepeatMasker on the reference
   6. CONSENSUS: Fasta file with CONSENSUS sequences for mobile elements in human
   7. REFERENCE: Fasta file for the reference human sequence
   8. SAMPLEID: Sample identified for naming the output VCF file

# Input files
* Bed files for VNTR, EXONS and REPEAT annotations can be downloaded from Zenodo for hg38 (link) and chm13 (link).
* Fasta containing consensus sequences for retroelements can be downloaded from Zenodo.
* TRF output can be generated as described bellow.

# Generation of TRF annotations for INS
1. Produce fasta with inserted sequences
python scripts/ins2fasta.py ins.vcf outDir

2. Execute TRF
trf insertions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs 1> ins_trf.out

# Generation of TRF annotations for DEL
1. Produce fasta with deleted sequences
python scripts/del2fasta.py del.vcf outDir

2. Execute TRF
trf deletions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs 1> del_trf.out

### Produce fasta with inserted sequences

## Execution
# Execution for INS (chm13):
python SVAN-INS.py ins.vcf ins_trf.out VNTR_chm13.bed EXONS_chm13.bed REPEATS_chm13.bed CONSENSUS.fa chm13.fa SAMPLEID

# Execution for DEL (chm13):
python SVAN-DEL.py del.vcf del_trf.out VNTR_chm13.bed EXONS_chm13.bed REPEATS_chm13.bed CONSENSUS.fa chm13.fa SAMPLEID 

## Output
SVAN produces as output a standard VCF file with SV annotations incorporated into the INFO field per each variant

## License
SVAN is distributed under GPL-3.0 License.

## Contact
Please open a case on the Github page for problems.

You may also contact Bernardo Rodriguez-Martin directly (e-mail omitted to reduce SPAM).
