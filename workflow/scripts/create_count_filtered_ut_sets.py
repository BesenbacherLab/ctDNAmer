import pandas as pd
import time
import os
import sys

## Helper functions ##
def get_gc_content(data):
    start = time.time()
    GC_content = []
    for index, row in data.iterrows():
        kmer = row["kmer"]
        GC = 0
        for base in kmer:
            if base in ["G", "C"]:
                GC += 1
        GC_perc = round(GC*100/len(kmer))
        GC_content.append(GC_perc)
    
    end = time.time()
    print(f'Time elapsed: {end - start}')
    sys.stdout.flush()
    return(GC_content)


# log and printing
os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
sys.stderr = open(snakemake.log[0], "a+")
sys.stdout = sys.stderr

# output
os.makedirs(os.path.dirname(snakemake.output.count_filtered_dir), exist_ok=True)
output_prefix = snakemake.params.out_pref 
print(f"Prefix for output files: {output_prefix}")

# input parameters
min_Tcount = int(snakemake.params.min_Tcount)
max_Tcount = int(snakemake.params.max_Tcount)
gc_lower = int(snakemake.params.gc_lower)
gc_upper = int(snakemake.params.gc_upper)
UT_size_upper_limit = int(snakemake.params.UT_size_upper_limit)
UT_min_size = int(snakemake.params.UT_min_size)

# input k-mer data
data = pd.read_csv(snakemake.input.UT_kmers, sep = '\t', header = None)
data.columns = ['kmer', 'count']
data['kmer'] = pd.Series(data['kmer'], dtype="string")
print(f'Data shape after reading it in: {data.shape}')

# filter based on the baseline minimum count
data_ci = data.loc[(data['count'] >= min_Tcount)]
data_ci.reset_index(drop = True, inplace = True)
print(f'Data shape after ci = {min_Tcount} filtering: {data_ci.shape}')

# filter based on the baseline maximum count
data_ci_cx = data_ci.loc[(data_ci['count'] <= max_Tcount)]
data_ci_cx.reset_index(drop = True, inplace = True)
print(f'Data shape after cx = {max_Tcount} filtering: {data_ci_cx.shape}')

# Calculate GC content
print(f'Calculating GC-content')
sys.stdout.flush()
GC_content = get_gc_content(data_ci_cx)
data_ci_cx["gc_content"] = GC_content

# filter based on GC content cutoffs
data_ci_cx_GC = data_ci_cx.loc[(data_ci_cx['gc_content'] >= gc_lower) & (data_ci_cx['gc_content'] <= gc_upper)]
print(f'Data shape after GC filtering: {data_ci_cx_GC.shape}')
print(f'Removed {data_ci_cx_GC.shape[0] - data_ci_cx_GC.shape[0]} with GC filtering')
data_ci_cx_GC.reset_index(drop = True, inplace = True)
sys.stdout.flush()

# ensure that the baseline UT set size is below the upper limit
cutoff = min_Tcount
while data_ci_cx_GC.shape[0] > UT_size_upper_limit:
    cutoff += 1
    print(f"Tumor count lower cutoff increased by one. Cutoff now is: {cutoff}")
    sys.stdout.flush()
    data_ci_cx_GC = data_ci_cx_GC.loc[(data_ci_cx_GC['count'] >= cutoff)]
data_ci_cx_GC.reset_index(drop = True, inplace = True)

# save baseline UT set
print("Finished baseline filtering.")
print(f'The baseline tumor count cutoff is: {cutoff}')
print(f'The baseline UT-kmers data set size is: {data_ci_cx_GC.shape[0]}')
os.makedirs(output_prefix, exist_ok=True)
data_ci_cx_GC.to_csv(f"{output_prefix}UT_kmers_filtered_minct{cutoff}.txt", sep="\t", header = False, index=False)
result_metadata = {'nUT_final': [data_ci_cx_GC.shape[0]],
                   'final_lower_cutoff': [cutoff], 
                   'final_upper_cutoff': [max_Tcount]}
result_metadata_df = pd.DataFrame.from_dict(result_metadata)
result_metadata_df.to_csv(f"{output_prefix}UT_mdata_minct{cutoff}.txt", sep="\t", header = True, index=False)

print("Starting count-based filtering")
sys.stdout.flush()

n_cap_datasets = 0
cutoffs_list = []
while data_ci_cx_GC.shape[0] > UT_min_size: # decrease UT set size by removing k-mers with the lowest count until the minimum size is reached
    n_cap_datasets += 1 # number of count-filtered data sets
    cutoff += 1 # current count cutoff
    cutoffs_list.append(cutoff)
    data_ci_cx_GC = data_ci_cx_GC.loc[(data_ci_cx_GC['count'] >= cutoff)] # filter data set based on count cutoff
    print(f"Tumor count lower cutoff increased by one for size reduction. The count ctuoff is now at : {cutoff}. The dataset size is now : {data_ci_cx_GC.shape[0]}")
    data_ci_cx_GC.reset_index(drop = True, inplace = True)
    sys.stdout.flush()
    
    data_ci_cx_GC.to_csv(f"{output_prefix}UT_kmers_filtered_minct{cutoff}.txt", sep="\t", header = False, index=False)
    
    # write count-filtered data set
    result_metadata = {'nUT_final': [data_ci_cx_GC.shape[0]],
                       'final_lower_cutoff': [cutoff], 
                       'final_upper_cutoff': [max_Tcount]}
    result_metadata_df = pd.DataFrame.from_dict(result_metadata)
    result_metadata_df.to_csv(f"{output_prefix}UT_mdata_minct{cutoff}.txt", sep="\t", header = True, index=False)


print(f"Number of capped datasets: {n_cap_datasets}")
print(f"Cutoff list: {cutoffs_list}")
    
print("Done!")
    
    
