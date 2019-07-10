import sys
import glob
import gzip
import subprocess
import string
import numpy

# This script will extract specific

# Usage: run ampliconseq_genotyper.py from directory containing trimmed sample directories. Results are printed to screen (stdout) 
# 1) update samples variable to select one file per directory. ampliconseq_genotyper.py should be run from the directory containing trimmed sample directories.
# 2) update left and right flanking sequences below to your specific locus of interest. The left flank must be directly flanking the genotyped sequence.
#   Note: longer sequences may result in higher specificity, but less sensitivity because imperfect matches are thrown out.
# 3) update base_pairs_genotyped to how many bps you want to examine in your data
# 4) update min_percent_coverage_per_allele. A lower percent (e.g. 0.1) will report more allelic variants and a higher percent reports only the most common ones.

#This selects one file from each sample directory as input (in this case pe1)
samples = glob.glob("*/trim_pe1*")

left_flank = "GATAGGAAAC"
right_flank = "GTCGTAAGT"

base_pairs_genotyped = 3 #I want the entire codon for kdr and ACE 
min_percent_coverage_per_allele = 0.1 # I only want to report the frequencies of the major alleles for readability.

# make reverse complement
old_chars = "ACGT"
replace_chars = "TGCA"
tab = string.maketrans(old_chars,replace_chars)

# This function extracts bps of interest
def get_genotype(trimmed_read,L_flanking_sequence,R_flanking_sequence,bps): # flanking sequences ensure read quality, but may also introduce some selection bias.
    RC_read = trimmed_read.translate(tab)[::-1] # check the reverse complement of all reads
    if L_flanking_sequence in trimmed_read and R_flanking_sequence in trimmed_read:
        genotype_pos = trimmed_read.find(L_flanking_sequence)
        genotype = trimmed_read[genotype_pos+len(L_flanking_sequence):genotype_pos+len(L_flanking_sequence)+bps]
    elif L_flanking_sequence in RC_read and R_flanking_sequence in RC_read:
        genotype_pos = RC_read.find(L_flanking_sequence)
        genotype = RC_read[genotype_pos+len(L_flanking_sequence):genotype_pos+len(L_flanking_sequence)+bps]
    else:
        genotype ="nada" # If the flanking sequences do not match, report "nada" (nothing), which is skipped later. nada reads are common.
    return genotype



# Output a list of the most common genotype sequences. L = left flanking sequencing. R = Sequence after SNP or codon
# This is important because I want to make a dataframe with the raw data but want to exclude super rare variants.
def make_genotype_list(glob_trimmed_read_file_per_sample, L, R, basepairs,allele_freq_cutoff):
    # make a dictionary with your data
    all_codons={}
    total_reads = 0
    reads_per_sample_list = []
    for sample in glob_trimmed_read_file_per_sample:
        sample_ID = sample.split('/')[0]
        all_sample_files = glob.glob("%s/trim_*" % (sample_ID)) # get all the trimmed data files. Run from directory containing trimmed sample directories.
        switch=0
        reads_per_sample=0
        for f in all_sample_files:
            for line in gzip.open(f):
                i=line.strip().split()
                if line[0]=="@":
                    switch=1
                    continue
                if switch:
                    switch=0 # I only want the next row so turn switch off
                    read=line.strip()
                    codon = get_genotype(read,L,R,basepairs)
                    if codon != "nada":
                        if codon not in all_codons:
                            all_codons[codon]=0
                        all_codons[codon]+=1
                        total_reads+=1
                        reads_per_sample+=1
        reads_per_sample_list.append(reads_per_sample)
    
    #calculate coverage threshold
    allele_freq_list = []
    gt_list = []
    for gt in all_codons:
        allele_freq_list.append(all_codons[gt])
        gt_list.append(gt)
    threshold = numpy.mean(allele_freq_list) * allele_freq_cutoff 
    GT_list = []
    for gt in all_codons:
        if all_codons[gt] > threshold:
            GT_list.append(gt)
    return [GT_list,numpy.mean(reads_per_sample_list)]


genotype_info = make_genotype_list(samples, left_flank, right_flank, base_pairs_genotyped, min_percent_coverage_per_allele)
codon_list = genotype_info[0]
ave_depth = genotype_info[1]

#print out a header:
print ",".join(['sample'] + codon_list)

for sample in samples:
    sample_ID = sample.split('/')[0]
    all_sample_files = glob.glob("%s/trim_*" % (sample_ID)) # get all the trimmed data files
    D={}
    switch=0
    total_reads=0
    for f in all_sample_files:
        for line in gzip.open(f):
            i=line.strip().split()
            if line[0]=="@":
                switch=1
                continue
            if switch:
                switch=0 # I only want the next row
                read=line.strip()

                codon = get_genotype(read,left_flank,right_flank,3)
                if codon != "nada":
                    if codon not in D:
                        D[codon]=0
                    D[codon]+=1
                    total_reads+=1
    
    if total_reads < ave_depth*0.1: # Here is a quality threshold. 10% of the mean sequencing depth.
        continue
    genotype_list=[sample_ID]
    codons=[sample_ID]
    for CODON in codon_list: # keep codon order for generation of a dataframe later
        codons.append(CODON)
        if CODON not in D:
            genotype_list.append('0')
            continue
        genotype_list.append(str(D[CODON]))
    # print out reads per allelic variant:
    print ",".join(genotype_list)
