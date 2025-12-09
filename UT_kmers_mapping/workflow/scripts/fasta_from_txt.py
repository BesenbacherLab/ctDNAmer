import os
import pandas as pd
import sys


# log and printing
os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
sys.stderr = open(snakemake.log[0], "a+")
sys.stdout = sys.stderr

# input file
kmers_f = snakemake.input.kmers

# output 
outfile = snakemake.output.fasta
os.makedirs(os.path.dirname(outfile), exist_ok=True)

# read in data
d = pd.read_csv(kmers_f, sep = "\t", header = None, names = ["kmer", "tumor", "gc"])
print("Head of the input kmers")
print(d.head())
print(d.shape[0])

with open(outfile, "w") as output:
    for index, row in d.iterrows():
        
        if index % 1000 == 0:
            print(index)
            sys.stdout.flush()
            
        header = f">kmer_{index}_tumor_{row["tumor"]}_gc_{row["gc"]}"
        output.write(header)
        output.write("\n")
        output.write(row["kmer"])
        output.write("\n")

print("Done!")    