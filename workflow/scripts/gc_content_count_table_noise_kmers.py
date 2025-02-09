import sys
import os
from collections import Counter

# logging, printing and output directory
os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.count_table), exist_ok=True)
sys.stderr = open(snakemake.log[0], "a+")
sys.stdout = sys.stderr

# input parameters
max_count = int(snakemake.params.max_count)
kmer_len = snakemake.params.k
print(f"Maximum count: {max_count}")
print(f"K-mer length: {kmer_len}")

# input data
file = open(snakemake.input.kmers, 'r')

counter = Counter()
index_c = 0
while True:
    line = file.readline()
    if not line:
        break
    
    index_c += 1
    if index_c % 100000 == 0: # track progress
        print(index_c)
        sys.stdout.flush()

    data = line.strip().split(" ")
    kmer = data[0]
    UT_count = int(data[1])
    donor_count = int(data[2])
    if UT_count <= max_count:
        GC = 0
        for i in kmer:
            if i in ["G", "C"]: # count GC bases
                GC += 1

        GC_perc = round(GC*100/kmer_len) # calculate GC content
        counter[(GC_perc, UT_count, donor_count)] += 1

file.close()

columns = ['GC','UT_count','donor_count', 'n_kmers']
outfile = open(snakemake.output.count_table, 'w')
sys.stdout = outfile # direct output to the correct file

# write output
print('\t'.join(columns))
for tup, count in counter.items():
    print('\t'.join(str(x) for x in tup) + '\t' + str(count))

sys.stdout = sys.stderr
outfile.close()
print("Done")
