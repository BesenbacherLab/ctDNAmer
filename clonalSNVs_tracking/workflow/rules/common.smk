import os
import sys
from pathlib import Path
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version


min_version("8.0.0")


###### Config file and sample sheets #####
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")


samples = pd.read_table(config["samples"]).set_index("sample_ID", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")


units = pd.read_table(config["samples_cfDNA"]).set_index("cfDNA_ID", drop=False)
units = units[units['sample_ID'].isin(samples.index.values.tolist())]
validate(units, schema="../schemas/units.schema.yaml")


###### Reference genome and VEP cahce info ######
reference_name = config["reference"]["name"]
reference_fasta = config["reference"]["fasta"]
vep_cache = config["vep_cache_path"]


###### Helper functions ######
def get_mutect_input(wildcards):
    """
    Mutect input. Read from sample sheet, return as a dict. 
    """
    input_dict = {"m_calls": samples.loc[wildcards.pt]["mutect"]}
    return input_dict


def get_sequenza_wiggle_track_from_ref_fasta():
    """
    Get the wiggle track path
    """
    basename_wo_suffix = ".".join(os.path.basename(reference_fasta).split(".")[:-1])
    output_name = f"results/sequenza/{basename_wo_suffix}.gc50Base.wig.gz"
    return output_name


def get_normal_bam(wildcards):
    """
    Germline BAM file
    """
    return samples.loc[wildcards.pt]["germline_bam"]


def get_tumor_bam(wildcards):
    """
    Tumor BAM file
    """
    return samples.loc[wildcards.pt]["tumor_bam"]


def get_cfDNA_bam(wildcards):
    """
    cfDNA BAM file
    """
    return units.loc[wildcards.cfDNA_ID]["cfDNA"]


def get_tracked_SNVs_input(wildcards):
    """
    Combine tracked clonal SNV data across cfDNA samples
    """
    SNVs_all_cf_pt = []
    units_pt = units[units['sample_ID'].isin([wildcards.pt])]
    units_cfDNA_id_list = units_pt.index.values.tolist()
    for unit in units_cfDNA_id_list: 
        SNVs_all_cf_pt.append(f"results/{wildcards.pt}/{unit}/SNVs_in_cfDNA.txt")
    return {"SNVs": SNVs_all_cf_pt}

def get_cfDNA_timepoints(wildcards):
    """
    Combine cfDNA samples' timepoints
    """
    SNVs_all_cf_pt = []
    units_pt = units[units['sample_ID'].isin([wildcards.pt])]
    units_pt_sub = units_pt[["timepoint"]]
    timepoints_dict = units_pt_sub.to_dict()["timepoint"]
    return timepoints_dict