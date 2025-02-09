import os
import pysam


# read in BAM file, search for index; if not found create index
def open_bam_create_index(bam_file):
    bai_filename1 = f"{bam_file}.bai"
    bai_filename2 = bam_file[:-1] + 'i'
    if not os.path.exists(bai_filename1) and not os.path.exists(bai_filename2):
        pysam.index(bam_file)
    return pysam.AlignmentFile(bam_file, "rb")


# Read class; defines query position information and read level stats
class Read: 
    def __init__(self, pileup_read):
        self.pos = pileup_read.query_position
        self.allel = pileup_read.alignment.query_sequence[self.pos]
        self.base_qual = pileup_read.alignment.query_qualities[self.pos]
        self.enddist = min(self.pos, len(pileup_read.alignment.query_sequence)-self.pos)
        
        self.secondary = pileup_read.alignment.is_secondary
        self.supplementary = pileup_read.alignment.is_supplementary
        self.mapq = pileup_read.alignment.mapping_quality
        
        cigar_stats = pileup_read.alignment.get_cigar_stats()[0]
        self.has_indel = sum(cigar_stats[1:4]) != 0
        self.has_clip = sum(cigar_stats[4:6]) != 0
        self.NM = cigar_stats[10]