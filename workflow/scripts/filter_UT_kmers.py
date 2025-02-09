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
os.makedirs(os.path.dirname(snakemake.output.filtered), exist_ok=True)
with open(snakemake.output.filtered, 'w') as fp:
    pass

# input parameters
gc_lower = int(snakemake.params.gc_lower)
gc_upper = int(snakemake.params.gc_upper)
max_Tcount = int(snakemake.params.max_Tcount)
min_data_length = int(snakemake.params.min_UT_size)

# input data
tumor_mean = pd.read_csv(snakemake.input.tumor_mean)
mean_count = tumor_mean["mean_count"][0].item()
round_tumor_mean = round(mean_count)
sd_count = tumor_mean["sd_count"][0].item()
print(f"Tumor mean, rounded: {round_tumor_mean}")
print(f"Tumor sd: {sd_count}")

# baseline count cutoffs based on tumor mean and sd
start_lcc = round(mean_count - sd_count)
start_ucc = round(mean_count + 5*sd_count)
print(f"Tumor starting cutoffs: {start_lcc}, {start_ucc}")

# input k-mer data
data = pd.read_csv(snakemake.input.UT_kmers, sep = '\t', header = None)
data.columns = ['kmer', 'count']
data['kmer'] = pd.Series(data['kmer'], dtype="string")
nUT_ci3_cx100 = data.shape[0]
print(f'Data shape after reading it in: {nUT_ci3_cx100}')

# calculate GC content
print(f'Calculating GC-content')
sys.stdout.flush()
GC_content = get_gc_content(data)
data["gc_content"] = GC_content

# write unfiltered UT set count summary with GC content information
data_write = data.drop(['kmer'], axis=1)
data_write = data_write.groupby(['count', 'gc_content'], as_index=False).size()
data_write.to_csv(snakemake.output.UT_data_gc_table, sep="\t", header = True, index=False)
sys.stdout.flush()

# filter based on GC content
data_GC = data.loc[(data['gc_content'] >= gc_lower) & (data['gc_content'] <= gc_upper)]
nUT_ci3_cx100_GC = data_GC.shape[0]
print(f'Data shape after GC filtering: {nUT_ci3_cx100_GC}')
print(f'Removed {nUT_ci3_cx100 - nUT_ci3_cx100_GC} with GC filtering')
sys.stdout.flush()
data_GC.reset_index(drop = True, inplace = True)

# filter based on baseline lower count cutoff
data_GC_ci = data_GC.loc[(data_GC['count'] >= start_lcc)]
data_GC_ci.reset_index(drop = True, inplace = True)
print(f'Data shape after lower count filtering: {data_GC_ci.shape}')
print(data_GC_ci.head())
sys.stdout.flush()

# filter based on baseline upper count cutoff
data_GC_ci_cx = data_GC_ci.loc[(data_GC_ci['count'] <= start_ucc)]
data_GC_ci_cx.reset_index(drop = True, inplace = True)
print(f'Data shape after upper count filtering: {data_GC_ci_cx.shape}')
print(data_GC_ci_cx.head())
sys.stdout.flush()

# UT set size after baseline count cutoff filtering
nUT_ci3_cx100_GC_BLcount = data_GC_ci_cx.shape[0]
print(f'Data shape after baseline count filtering: {nUT_ci3_cx100_GC_BLcount}')

final_lcc = start_lcc
final_ucc = start_ucc
lower_cutoff_limit = 2

# filtering path trackers
above_min = False
increasing_lcc = False
decreasing_ucc = False
below_minimum = False

# UT set size above target 
if data_GC_ci_cx.shape[0] > min_data_length:
    above_min = True
    while data_GC_ci_cx.shape[0] > min_data_length and final_lcc < (round_tumor_mean-1): # increase lower cutoff
        increasing_lcc = True
        final_lcc += 1
        data_GC_ci_cx = data_GC_ci_cx.loc[(data_GC_ci_cx['count'] >= final_lcc)]
        
        print(f"Tumor count lower cutoff increased by one. Cutoff now is: {final_lcc}")
        sys.stdout.flush()

    if data_GC_ci_cx.shape[0] > min_data_length and final_lcc == (round_tumor_mean-1): 
        while data_GC_ci_cx.shape[0] > min_data_length and final_ucc > (round_tumor_mean+1): # decrease upper cutoff, if still above target size
            increasing_lcc = False
            decreasing_ucc = True
            
            final_ucc -= 1
            data_GC_ci_cx = data_GC_ci_cx.loc[(data_GC_ci_cx['count'] <= final_ucc)]

            print(f"Tumor count upper cutoff decreased by one. Cutoff now is: {final_ucc}")
            sys.stdout.flush()

# UT set size below target 
elif data_GC_ci_cx.shape[0] < min_data_length:
    while data_GC_ci_cx.shape[0] < min_data_length and final_ucc < 100: # increase upper cutoff

        final_ucc += 1
        print(f"Tumor count upper cutoff increased by one. Cutoff now is: {final_ucc}")
        sys.stdout.flush()

        nrow_before = data_GC_ci_cx.shape[0]
        higher_cutoff_kmers = data_GC.loc[(data_GC['count'] == final_ucc)]
        data_GC_ci_cx = pd.concat([data_GC_ci_cx, higher_cutoff_kmers])

        print(f'Data shape after increasing upper cutoff: {data_GC_ci_cx.shape}')
        print(f'Added {data_GC_ci_cx.shape[0] - nrow_before} with upper cutoff increasing')
        sys.stdout.flush()

    if data_GC_ci_cx.shape[0] < min_data_length and final_ucc == 100:
        while data_GC_ci_cx.shape[0] < min_data_length and final_lcc > 3: # decrease lower cutoff, if still below target size
            final_lcc -= 1
            print(f"Tumor count lower cutoff decreased by one. Cutoff now is: {final_lcc}")
            sys.stdout.flush()

            nrow_before = data_GC_ci_cx.shape[0]
            lower_cutoff_kmers = data_GC.loc[(data_GC['count'] == final_lcc)]
            data_GC_ci_cx = pd.concat([data_GC_ci_cx, lower_cutoff_kmers])
    
            print(f'Data shape after decreasing lower cutoff: {data_GC_ci_cx.shape}')
            print(f'Added {data_GC_ci_cx.shape[0] - nrow_before} with lower cutoff decreasing')
            sys.stdout.flush()

shifting_3count_set = False
final_size_adjustment_w_ucc = False

if above_min and data_GC_ci_cx.shape[0] < min_data_length: # if we decreased the UT set size before below the target, adjust it back to above target
    print(f"UT set size now below the target size {min_data_length}, current cutoffs: {final_lcc}, {final_ucc}")
    print("Cutoffs will be now adjusted by one, so the final data set size is (just) above the target size.")
    if increasing_lcc: # adjusting the lower cutoff back by one, if this was the cutoff last adjusted above
        final_lcc -= 1
        print(f"Tumor count lower cutoff decreased by one. Cutoff now is: {final_lcc}")

        nrow_before = data_GC_ci_cx.shape[0]
        lower_cutoff_kmers = data_GC.loc[(data_GC['count'] == final_lcc)]
        data_GC_ci_cx = pd.concat([data_GC_ci_cx, lower_cutoff_kmers])

        print(f'Data shape after decreasing lower cutoff: {data_GC_ci_cx.shape}')
        print(f'Added {data_GC_ci_cx.shape[0] - nrow_before} with lower cutoff decreasing')
        sys.stdout.flush()

        # if only lower bound increased before, checking now if upper bound can be decreased while still staying above the target size
        while data_GC_ci_cx.shape[0] > min_data_length and final_ucc > (round_tumor_mean+1): 
            final_size_adjustment_w_ucc = True
            final_ucc -= 1
            
            data_GC_ci_cx = data_GC_ci_cx.loc[(data_GC_ci_cx['count'] <= final_ucc)]
            print(f"Final size adjustment with ucc: upper cutoff decreased by one. Cutoff now is: {final_ucc}")

        if final_size_adjustment_w_ucc and data_GC_ci_cx.shape[0] < min_data_length: # if we decreased the UT set size before below the target, adjust it back to above target
            final_ucc += 1
            print(f"Final size adjustment with ucc, last step: upper cutoff increased (back) by one. Cutoff now is: {final_ucc}")
            nrow_before = data_GC_ci_cx.shape[0]
            higher_cutoff_kmers = data_GC.loc[(data_GC['count'] == final_ucc)]
            data_GC_ci_cx = pd.concat([data_GC_ci_cx, higher_cutoff_kmers])
    
            print(f'Data shape after increasing upper cutoff: {data_GC_ci_cx.shape}')
            print(f'Added {data_GC_ci_cx.shape[0] - nrow_before} with upper cutoff increasing')
            sys.stdout.flush()
    
    if decreasing_ucc: # adjusting the upper cutoff back by one, if this was the cutoff last adjusted above
        final_ucc += 1
        print(f"Upper cutoff increased (back) by one. Cutoff now is: {final_ucc}")
        
        nrow_before = data_GC_ci_cx.shape[0]
        higher_cutoff_kmers = data_GC.loc[(data_GC['count'] == final_ucc)]
        data_GC_ci_cx = pd.concat([data_GC_ci_cx, higher_cutoff_kmers])

        print(f'Data shape after decreasing lower cutoff: {data_GC_ci_cx.shape}')
        print(f'Added {data_GC_ci_cx.shape[0] - nrow_before} with upper cutoff increasing')
        sys.stdout.flush()

# if the data size (still) above target, check if increasing both the lower and upper cutoff by one still keeps it above target; keep increasing until boundary
if above_min and data_GC_ci_cx.shape[0] > min_data_length: 
    shifting_3count_set = True
    print(f'lower and upper cutoffs brought to limits, but size still above min target size. Increasing the lower and upper cutoffs together by one to decrease the data set size. Current cutoffs: {final_lcc}, {final_ucc}')
    while data_GC_ci_cx.shape[0] > min_data_length:
        final_lcc += 1
        final_ucc += 1
        data_GC_ci_cx = data_GC[data_GC['count'].between(final_lcc, final_ucc)]
        print(f'Increased both cutoffs by one, data set size now: {data_GC_ci_cx.shape[0]}')
        sys.stdout.flush()

# adjust back by one to ensure UT set size above target
if shifting_3count_set:
    final_lcc -= 1
    final_ucc -= 1
    data_GC_ci_cx = data_GC[data_GC['count'].between(final_lcc, final_ucc)]
    print(f'Decreased both cutoffs (back) by one, to get the final data set size: {data_GC_ci_cx.shape[0]}')
    sys.stdout.flush()

data_GC_ci_cx.reset_index(drop = True, inplace = True)

print("Finished tumor count cutoff filtering.")
print(f'The final tumor count cutoffs are: {final_lcc}, {final_ucc}')
print(f'The final UT-kmers data set shape is: {data_GC_ci_cx.shape[0]}')
sys.stdout.flush()
print(data_GC_ci_cx.head())

# write output
data_GC_ci_cx.to_csv(snakemake.output.filtered, sep="\t", header = False, index=False)

result_metadata = {'nUT_ci3_cx100': [nUT_ci3_cx100],
                   'nUT_ci3_cx100_GC': [nUT_ci3_cx100_GC],
                   'nUT_ci3_cx100_GC_BLcount': [nUT_ci3_cx100_GC_BLcount],
                   'nUT_final': [data_GC_ci_cx.shape[0]],
                   'start_lower_cutoff': [start_lcc],
                   'start_upper_cutoff': [start_ucc], 
                   'final_lower_cutoff': [final_lcc], 
                   'final_upper_cutoff': [final_ucc], 
                   'tumor_mean_count': [mean_count],
                   'rounded_tumor_mean_count': [round_tumor_mean]}
result_metadata_df = pd.DataFrame.from_dict(result_metadata)
result_metadata_df.to_csv(snakemake.output.filtered_mdata, sep="\t", header = True, index=False)

print("Done!")
