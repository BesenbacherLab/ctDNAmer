import os
import pandas as pd
import sys
import pysam
import statistics
from pysam_utils import open_bam_create_index, Read

# log and printing
os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
sys.stderr = open(snakemake.log[0], "a+")
sys.stdout = sys.stderr

# input files and parameters
bam_file = snakemake.input.bam_file
snvs_file = snakemake.input.snvs_to_track
fasta_file = snakemake.input.reference_fasta

# output 
outfile = snakemake.output.outfile
os.makedirs(os.path.dirname(outfile), exist_ok=True)

# read in data
bam = open_bam_create_index(bam_file)
snvs = pd.read_csv(snvs_file, sep = "\t")
print("Head of the input snvs")
print(snvs.head())
print(snvs.shape[0])

print_prog = False
n_snvs_pass = 0

with open(outfile, "w") as output:
    output.write("chrom\tstart\tend\tref\talt\ttumor_DP\ttumor_n_ref\ttumor_n_alt\tcfDNA_DP_total\tcfDNA_n_ref_total\tcfDNA_n_alt_total\tcfDNA_DP_filt\tcfDNA_n_ref_filt\tcfDNA_n_alt_filt\n")
    for index, row in snvs.iterrows(): # iterate over input SNVs
        if index % 100 == 0: # progress tracker
            print_prog = True
        else: 
            print_prog = False
        
        # save SNV information
        chr_i = row["chrom"]
        start = row["start"]
        end = row["end"]
        ref_i = row["ref"]
        ref_i = ref_i.upper()
        alt_i = row["alt"]
        alt_i = alt_i.upper()
        tumor_DP = row["DP"]
        tumor_n_ref = row["n_ref"]
        tumor_n_alt = row["n_alt"]

        for pileupcolumn in bam.pileup(str(chr_i), int(start), int(end), truncate = True, ignore_overlaps = True, min_base_quality = 0, compute_baq = False):
            if print_prog:
                print(f'SNP nr {index} --> Chr: {chr_i}, from: {start}, SNV: {ref_i}->{alt_i}')
                sys.stdout.flush()
            position = [x.upper() for x in pileupcolumn.get_query_sequences(add_indels = True)] # all bases in the current position
            alt_bases = set([x for x in position if x in ['A','C','G','T'] if x != ref_i]) # set of altenrative alleles in the current position
            if len(alt_bases) > 1: # ignore variant if more than 1 unique alternative allele in cfDNA
                continue
            if len(alt_bases) == 1: 
                alt = list(alt_bases)[0]
                if alt.upper() != alt_i.upper():
                    continue

            n = 0
            n_ref_total = 0
            n_alt_total = 0
            n_high_qual_ref = 0
            n_high_qual_alt = 0
            n_high_qual = 0
            basequals_list = []
            enddist_list = []
            for pileup_read in pileupcolumn.pileups:
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue
                read = Read(pileup_read)
                n += 1
                # read level quality filter: 0 indels and clips; only 1 ref-read mismatch across read; not a secondary or supplementary alingment; MQ = 60
                if read.has_indel or \
                    read.has_clip or \
                    read.NM > 1 or \
                    read.secondary == True or \
                    read.supplementary == True or \
                    read.mapq < 60:
                    if read.allel.upper() not in [ref_i, alt_i]:
                        continue
                    else:
                        if read.allel.upper() == ref_i:
                            n_ref_total +=1
                        elif read.allel.upper() == alt_i:
                            n_alt_total += 1
                else:
                    basequals_list.append(read.base_qual)
                    enddist_list.append(read.enddist)
                    n_high_qual += 1 # if all above conditions met, read is classified as high-quality
                    if read.allel.upper() == ref_i:
                        n_ref_total +=1
                        n_high_qual_ref +=1
                    elif read.allel.upper() == alt_i:
                        n_alt_total += 1
                        n_high_qual_alt += 1
            
            # position level quality filter: 90% reads high quality, high-quality coverage >= 15, median BQ > 20, median distance to read end > 5
            if n != 0: 
                if n_high_qual/n >= 0.9 and \
                    n_high_qual >= 15 and \
                    statistics.median(basequals_list) > 20 and \
                    statistics.median(enddist_list) > 5:
                    # write quality filtered output
                    n_snvs_pass += 1
                    output.write(f"{chr_i}\t{start}\t{end}\t{ref_i}\t{alt_i}\t{tumor_DP}\t{tumor_n_ref}\t{tumor_n_alt}\t{n}\t{n_ref_total}\t{n_alt_total}\t{n_high_qual}\t{n_high_qual_ref}\t{n_high_qual_alt}\n") 

print(f"Out of {snvs.shape[0]} clonal SNVs, {n_snvs_pass} passed in cfDNA")
print("Done!")    