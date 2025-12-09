import os
import pandas as pd
import sys
import pysam


# log and printing
os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
sys.stderr = open(snakemake.log[0], "a+")
sys.stdout = sys.stderr

# input file
bam_f = snakemake.input.kmers
mutect_f = snakemake.input.mutect
delly_f = snakemake.input.delly

# output 
outfile = snakemake.output.summary
os.makedirs(os.path.dirname(outfile), exist_ok=True)

one_mm_outfile = snakemake.output.one_mm_out
os.makedirs(os.path.dirname(one_mm_outfile), exist_ok=True)

more_unique_outfile = snakemake.output.more_unique_out
os.makedirs(os.path.dirname(more_unique_outfile), exist_ok=True)


class Read:
    def __init__(self, pileup_read):
        
        self.ref = pileup_read.reference_name
        self.start = pileup_read.reference_start
        self.end = pileup_read.reference_end

        self.query_name = pileup_read.query_name
        self.query_sequence = pileup_read.query_sequence
        self.is_unmapped = pileup_read.is_unmapped

        self.mapq = pileup_read.mapping_quality
        
        cigar_stats = pileup_read.get_cigar_stats()[0]
        self.has_indel = sum(cigar_stats[1:4]) != 0
        self.has_clip = sum(cigar_stats[4:6]) != 0
        self.NM = cigar_stats[10]

# read in data
bam = pysam.AlignmentFile(bam_f,'rb')
mutect = pysam.VariantFile(mutect_f)
delly = pysam.VariantFile(delly_f)


i = 0
with open(outfile, "w") as output, open(one_mm_outfile, "w") as one_mm_output, open(more_unique_outfile, "w") as more_unique_output:
    output.write("kmer\tchr\talign_start\talign_end\tMQ\tis_unmapped\thas_indel\thas_clip\tn_mismatches\tn_mutect_matches\tmutect_variant_list\tn_delly_matches\tdelly_variant_list\n")
    for read in bam.fetch(until_eof=True):
        
        # fetch read information
        read_ob = Read(read)

        # kmer, tumor, gc
        kmer = read_ob.query_sequence
        kmer_info = read_ob.query_name
        tumor_count = kmer_info.split("_")[3]
        gc = kmer_info.split("_")[5]

        if (not read_ob.is_unmapped) & (read_ob.NM == 1) & (not read_ob.has_clip) & (not read_ob.has_indel): 
            one_mm_output.write(f"{kmer}\t{tumor_count}\t{gc}\n")
        else: 
            more_unique_output.write(f"{kmer}\t{tumor_count}\t{gc}\n")
            
        i += 1
        if i % 1000 == 0:
            print(i)
    
        n_mutect_overlap = 0
        mut_var_list = []

        n_delly_overlap = 0
        delly_var_list = []
        if not read_ob.is_unmapped:
            for var in mutect.fetch(read_ob.ref, start=read_ob.start, stop=read_ob.end):
                n_mutect_overlap += 1
                mut_var_str = f"{var.chrom}_{var.pos}_{var.ref}_{"/".join(list(var.alts))}_{",".join(var.filter.keys())}"
                mut_var_list.append(mut_var_str)

            for var in delly.fetch(read_ob.ref, start=read_ob.start, stop=read_ob.end):
                n_delly_overlap += 1
                del_var_str = f"{var.chrom}_{var.pos}_{var.ref}_{"/".join(list(var.alts))}_{",".join(var.filter.keys())}"
                delly_var_list.append(del_var_str)
                    
        output.write(f"{read_ob.query_name}\t{read_ob.ref}\t{read_ob.start}\t{read_ob.end}\t{read_ob.mapq}\t{read_ob.is_unmapped}\t{read_ob.has_indel}\t{read_ob.has_clip}\t{read_ob.NM}\t{n_mutect_overlap}\t{";".join(mut_var_list)}\t{n_delly_overlap}\t{";".join(delly_var_list)}\n")


delly.close()
mutect.close()
bam.close()

print("Done!")    