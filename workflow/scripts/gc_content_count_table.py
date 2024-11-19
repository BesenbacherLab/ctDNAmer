import sys
import os
from collections import Counter

# logging
os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.count_table), exist_ok=True)
sys.stderr = open(snakemake.log[0], "a+")
sys.stdout = sys.stderr


# input
max_count = int(snakemake.params.max_count)
kmer_len = int(snakemake.params.k)
print(f"Maximum count: {max_count}")
print(f"K-mer length: {kmer_len}")

counter = Counter()
index_c = 0
file = open(snakemake.input.kmers, 'r')

while True:
    line = file.readline()
    if not line:
        break
    
    index_c += 1
    if index_c % 100000 == 0:
        print(index_c)
        sys.stdout.flush()

    data = line.strip().split("\t")
    kmer = data[0]
    count = int(data[1])
    if count <= int(max_count):
        GC = 0
        for i in kmer:
            if i in ["G", "C"]:
                GC += 1

        GC_perc = round(GC*100/kmer_len)
        counter[(GC_perc, count)] += 1

file.close()

columns = ['GC','count', 'n_kmers']
outfile = open(snakemake.output.count_table, 'w')
sys.stdout = outfile
print('\t'.join(columns))
for tup, count in counter.items():
    print('\t'.join(str(x) for x in tup) + '\t' + str(count))

sys.stdout = sys.stderr
outfile.close()
print("Done")