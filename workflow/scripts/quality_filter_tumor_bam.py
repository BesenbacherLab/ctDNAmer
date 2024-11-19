import pysam
import sys
import pandas as pd
import os
    
    
def FilterReadsBam(in_file, out_file):
    print("Starting to filter")
    sys.stdout.flush()

    max_quality = 37
    bam_in = pysam.Samfile(in_file, 'rb')
    bam_out = pysam.Samfile(out_file, 'wb', template=bam_in)

    reads = 0
    passing_reads = 0
    filtered = 0
    filtered_out = 0
    for read in bam_in.fetch():
        reads += 1
        if reads % 10000000 == 0:
            print(reads)
            sys.stdout.flush()

        passing_read = True
        if len(read.seq) == len(read.qual):
            len_read = len(read.seq)
            low_qual_list = [i for i,x in enumerate(read.qual) if ord(x)-33 < max_quality]
            if len(low_qual_list) > 0:
                passing_read = False
            
            new_read = "".join(["N" if i in low_qual_list or i < 5 or i > (len_read - 6) else x for i, x in enumerate(read.seq)])
    
            read.seq = new_read
            if passing_read: 
                passing_reads += 1
            else:
                filtered += 1
    
            bam_out.write(read)
        else:
            filtered_out += 1


    print(f"Total reads analyzed: {reads}")
    print(f"Reads passing the filtering:  {passing_reads}")
    print(f"Reads where at least one base set to N:  {filtered}")
    print(f"Reads which were removed, as quality and sequence were of unequal length:  {filtered_out}")
    
    bam_in.close()
    bam_out.close()

    return "Finished"

# log and printing
os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
sys.stderr = open(snakemake.log[0], "a+")
sys.stdout = sys.stderr

# output
os.makedirs(os.path.dirname(snakemake.output.filtered_files), exist_ok=True)
with open(snakemake.output.filtered_files, 'w') as fp:
    pass

# input
in_files = snakemake.input.input_files
out_files = snakemake.output.filtered_files
output_path_prefix = os.path.dirname(out_files) + "/"


input_files = open(in_files, 'r')
count = 0
for line in input_files:
    count += 1
    input_filepath_i = line.strip()
    input_filename = os.path.basename(input_filepath_i)
    print(input_filename)
    file_prefix_i = '.'.join(input_filename.split(".")[:-1])
    print(file_prefix_i)
    
    out_list = []
    outfile_i = file_prefix_i + ".maxqual_filtered.bam"
    output_filepath_i = output_path_prefix + outfile_i
    status = FilterReadsBam(input_filepath_i, output_filepath_i)
    print(status + f" with file number {count}, ({input_filename})")
    out_list.append(output_filepath_i)
    sys.stdout.flush()


print("Writing filtered files list")
with open(out_files, 'w') as f:
    for line in out_list:
        f.write(f"{line}\n")

print("Done!")
    

