import pandas as pd
from pathlib import Path
import os
import sys
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("8.0.0")


###### Config file and sample sheets #####
configfile: "config/config_count_filter_test.yaml"
validate(config, schema="../../schemas/config_count_filter_test.schema.yaml")

# patient data
samples = pd.read_table(config["samples"]).set_index("sample_ID", drop=False)
validate(samples, schema="../../schemas/samples_count_filter_test.schema.yaml")

# cfDNA sample data
units = pd.read_table(config["samples_cfDNA"]).set_index("cfDNA_ID", drop=False)
validate(units, schema="../../schemas/units_count_filter_test.schema.yaml")

preop_units = units.loc[units["timepoint"] <= 0, :].copy()
postop_index = units.groupby(['sample_ID'])['timepoint'].transform("max") == units['timepoint']
postop_units = units[postop_index]
postop_units.set_index("sample_ID", drop=False)

# unmatched cfDNA data for empirical noise estimation
donors = pd.read_table(config["donors"]).set_index("donor_ID", drop=False)
validate(donors, schema="../../schemas/donors.schema.yaml")

# germline union samples
samples_glu = pd.read_table(config["samples_glu"]).set_index("sample_ID", drop=False)
validate(samples_glu, schema="../../schemas/samples_glu.schema.yaml")

glu_sample_nr = [str(x) for x in range(1, len(samples_glu.index.tolist()) + 1)]
glu_combinations = ["1_" + str(x) for x in range(2, len(samples_glu.index.tolist()) + 1)]
glu_codes_dict = dict(zip(glu_sample_nr, samples_glu.index.tolist()))

###### Helper functions #####

def k():
    """
    K-mer length
    """
    return config["k"]


def get_input(path_str, sample_list):
    """
    Handles input data files. 
    Writes a text file with input file names and retuns the text file as a dict.
    """
    res_path = f"results/{path_str}"
    Path(res_path).mkdir(parents=True, exist_ok=True)
    if not os.path.isfile(f"{res_path}input_files.txt"):
        with open(f"{res_path}input_files.txt", 'w') as file:
            file.write('\n'.join(sample_list))
    
    tmp_f = f"kmc_bins/count_kmers/{path_str}"
    Path(tmp_f).mkdir(parents=True, exist_ok=True)
    
    return {"input_files": f"{res_path}input_files.txt"}


def get_fileformat(file_list, sample_type, sample_ID):
    """
    Input data format for k-mer counting command.
    Either BAM ("BAM") or FASTQ ("FASTQ") format allowed
    """
    file_suf = set()
    for f in file_list:
        file_suf.add(f.split(".")[-1])
    if len(file_suf) > 1:
        sys.exit(f'More than 1 type of input files for {sample_type} data; sample_ID: {sample_ID}')
    else:
        file_suffix = list(file_suf)[0]
        if file_suffix.lower() == "bam":
            file_format = "bam"
        else:
            file_format = "q"
        return file_format


def get_germline_union_input(wildcards):
    """
    Germline union input files
    """
    path_str = f"germline_union/{wildcards.glu_type}/{wildcards.glu_sample}/"
    sample_list = samples_glu.loc[wildcards.glu_sample]["data_files"].split(",")
    return get_input(path_str, sample_list)


def get_glu_fileformat(wildcards):
    """
    Germline union input format
    """
    file_list = samples_glu.loc[wildcards.glu_sample]["data_files"].split(",")
    return get_fileformat(file_list, "germline union", wildcards.glu_sample)


def get_glu_input(wildcards):
    """
    Germline union creation order up to the last combination
    """
    input_file_prefix_list = recurse_glu_samples(wildcards)
    l_pre = [file_prefix + ".kmc_pre" for file_prefix in input_file_prefix_list] 
    return l_pre

def get_glu_input_last(wildcards):
    """
    Germline union creation last combination
    """
    input_file_prefix_list = get_final_glu_prefix(wildcards)
    l_pre = [file_prefix + ".kmc_pre" for file_prefix in input_file_prefix_list]
    return  l_pre

def get_final_glu_prefix(wildcards):
    """
    Germline union creation last combination files
    """
    new_sample = glu_codes_dict[str(wildcards.glu_samples_comb_last.split("_")[1])]
    dtype_new = samples_glu.loc[new_sample]["data_type"]
    prev_union = "1_" + str(int(wildcards.glu_samples_comb_last.split("_")[1])-1)

    return [f"results/germline_union/germline_union_{prev_union}", 
            f"results/germline_union/{dtype_new}/{new_sample}/kmers"]

def recurse_glu_samples(wildcards):
    """
    Germline union creation. Recursion over input samples
    """
    if wildcards.glu_samples_comb == "1_2":
        sample1 = glu_codes_dict["1"]
        sample2 = glu_codes_dict["2"]
        dtype1 = samples_glu.loc[sample1]["data_type"]
        dtype2 = samples_glu.loc[sample2]["data_type"]
        return [f"results/germline_union/{dtype1}/{sample1}/kmers", f"results/germline_union/{dtype2}/{sample2}/kmers"]
    else:
        new_sample = glu_codes_dict[str(wildcards.glu_samples_comb.split("_")[1])]
        dtype_new = samples_glu.loc[new_sample]["data_type"]
        prev_union = "1_" + str(int(wildcards.glu_samples_comb.split("_")[1])-1)

        return [f"results/germline_union/germline_union_{prev_union}", 
                f"results/germline_union/{dtype_new}/{new_sample}/kmers"]


def get_final_glu():
    """
    Germline union path
    """
    return f"results/germline_union/final_germline_union_{glu_combinations[-1]}.kmc_pre"


def get_germline_input(wildcards):
    """
    Germline data input files
    """
    path_str = f"patients/{wildcards.pt}/germline/"
    sample_list = samples.loc[wildcards.pt]["germline"].split(",")
    return get_input(path_str, sample_list)


def get_germline_fileformat(wildcards):
    """
    Germline data input format
    """
    file_list = samples.loc[wildcards.pt]["germline"].split(",")
    return get_fileformat(file_list, "germline", wildcards.pt)
    

def get_tumor_input(wildcards):
    """
    Primary tumor data input files
    """
    path_str = f"patients/{wildcards.pt}/tumor/"
    sample_list = samples.loc[wildcards.pt]["tumor"].split(",")
    return get_input(path_str, sample_list)


def get_tumor_fileformat(wildcards):
    """
    Primary tumor data input format
    """
    file_list = samples.loc[wildcards.pt]["tumor"].split(",")
    return get_fileformat(file_list, "tumor", wildcards.pt)


def get_cfDNA_input(wildcards):
    """
    cfDNA data input files
    """
    path_str = f"patients/{wildcards.pt}/{wildcards.cfDNA_ID}/"
    sample_list = units.loc[wildcards.cfDNA_ID]["cfDNA"].split(",")
    return get_input(path_str, sample_list)


def get_cfDNA_fileformat(wildcards):
    """
    cfDNA data input format
    """
    file_list = units.loc[wildcards.cfDNA_ID]["cfDNA"].split(",")
    return get_fileformat(file_list, "cfDNA", wildcards.cfDNA_ID)


def get_donor_input(wildcards):
    """
    Unmatched cfDNA data input files
    """
    path_str = f"donors/{wildcards.donor}/"
    sample_list = donors.loc[wildcards.donor]["cfDNA"].split(",")
    return get_input(path_str, sample_list)


def get_donor_fileformat(wildcards):
    """
    Unmatched cfDNA data input format
    """
    file_list = donors.loc[wildcards.donor]["cfDNA"].split(",")
    return get_fileformat(file_list, "donor", wildcards.donor)

def aggregate_count_filtering_test_preop_estimates(wildcards):
    """
    Pre-treatment/baseline cfDNA sample model estimates
    """
    preop_sample = preop_units.index[preop_units.sample_ID == wildcards.pt].item()
    preop_estimates = f"results/patients/{wildcards.pt}/{preop_sample}/unique_tumor_count_filter_test/tf_estimation/minct{wildcards.ct_cutoff}/preop_tf_estimate.csv"
    return preop_estimates


def count_filtered_preop_estimates(wildcards):
    """
    Aggregation of pre-treatment/baseline TF estimates across count-filtered UT sets
    """
    checkpoint_out = checkpoints.create_count_filtered_ut_sets.get(**wildcards).output[0]

    return expand(
        "results/patients/{pt}/{preop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/preop_tf_estimate.csv",
        pt = wildcards.pt,
        preop_cfDNA_ID = preop_units.index[preop_units.sample_ID == wildcards.pt].item(),
        ct_cutoff = glob_wildcards(os.path.join(checkpoint_out, "UT_kmers_filtered_minct{ct_cutoff}.txt")).ct_cutoff,
    )


def count_filtered_postop_estimates(wildcards):
    """
    Aggregation of post-treatment/ctDNA-negative TF estimates across count-filtered UT sets
    """
    checkpoint_out = checkpoints.create_count_filtered_ut_sets.get(**wildcards).output[0] 

    return expand(
        "results/patients/{pt}/{postop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/postop_tf_estimate.csv",
        pt = wildcards.pt,
        postop_cfDNA_ID = postop_units.index[postop_units.sample_ID == wildcards.pt].item(),
        ct_cutoff = glob_wildcards(os.path.join(checkpoint_out, "UT_kmers_filtered_minct{ct_cutoff}.txt")).ct_cutoff,
    )
