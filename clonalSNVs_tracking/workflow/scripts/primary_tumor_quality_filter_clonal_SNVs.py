import os
import pandas as pd
import pysam
import sys
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

print_prog = False
n_snvs_pass = 0

with open(outfile, "w") as output:
    output.write("chrom\tstart\tend\tref\talt\tDP\tn_ref\tn_alt\n")
    for index, row in snvs.iterrows(): # iterating over each SNV
        if index % 100 == 0: # progress tracking
            print_prog = True
        else: 
            print_prog = False
        chr_i = row["chr"]
        start = row["from"]
        end = row["to"]
        ref_i = row["ref"]
        ref_i = ref_i.upper()
        alt_i = row["alt"]
        alt_i = alt_i.upper()
        for pileupcolumn in bam.pileup(str(chr_i), int(start), int(end), truncate = True, ignore_overlaps = True, min_base_quality = 0, compute_baq = False):
            if print_prog:
                print(f'SNP nr {index} --> Chr: {chr_i}, from: {start}, SNV: {ref_i}->{alt_i}')
                sys.stdout.flush()
    
            position = [x.upper() for x in pileupcolumn.get_query_sequences(add_indels = True)] # all bases in the current position
            alt_bases = set([x for x in position if x in ['A','C','G','T'] if x != ref_i]) # set of alternatives in current position
            if len(alt_bases) == 0 or len(alt_bases) > 1: # only 1 unique alternative allowed
                continue
            alt = list(alt_bases)[0] # alternative observed
            if alt.upper() != alt_i:
                continue
    
            n = 0
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
                    read.mapq < 60 or \
                    read.allel.upper() not in [ref_i, alt_i]:
                    continue
                else:
                    basequals_list.append(read.base_qual)
                    enddist_list.append(read.enddist)
                    n_high_qual += 1 # if all above conditions met, read is classified as high-quality
                    if read.allel.upper() == ref_i:
                        n_high_qual_ref +=1
                    elif read.allel.upper() == alt:
                        n_high_qual_alt += 1

            # position level quality filter: 90% reads high quality, high-quality coverage >= 20, median BQ > 20, median distance to read end > 5, 2 high quality alternative bases
            if n_high_qual/n >= 0.9 and \
                n_high_qual >= 20 and \
                statistics.median(basequals_list) > 20 and \
                statistics.median(enddist_list) > 5 and \
                n_high_qual_alt > 2:
                # write quality filtered output
                n_snvs_pass += 1
                output.write(f"{chr_i}\t{start}\t{end}\t{ref_i}\t{alt}\t{n_high_qual}\t{n_high_qual_ref}\t{n_high_qual_alt}\n")

print(f"Out of {snvs.shape[0]} clonal SNVs, {n_snvs_pass} passed for tracking in cfDNA")
print("Done!")