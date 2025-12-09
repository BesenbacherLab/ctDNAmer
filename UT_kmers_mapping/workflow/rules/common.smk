import os
import sys
from pathlib import Path
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version


min_version("8.0.0")


###### Config file and sample sheets #####
configfile: "config/config_UT_mapping_kmers.yaml"
#validate(config, schema="../schemas/config.schema.yaml")

# patient data
samples = pd.read_table(config["samples"]).set_index("sample_ID", drop=False)
#validate(samples, schema="../schemas/samples.schema.yaml")

# cfDNA data
units = pd.read_table(config["samples_cfDNA"]).set_index("cfDNA_ID", drop=False)
validate(units, schema="../schemas/units.schema.yaml")

preop_units = units.loc[units["timepoint"] <= 0, :].copy()
preop_units["preop_tf_path"] = preop_units.apply(lambda row: f"results/UT_kmers_mapping/patients/{row.sample_ID}/{UT_subset}/{row.name}/tf_estimation/preop_tf_estimate.csv", axis = 1)
postop_units = units.loc[units["timepoint"] > 0, :].copy()

# unmatched cfDNA data for empirical noise estimation
donors = pd.read_table(config["donors"]).set_index("donor_ID", drop=False)
#validate(donors, schema="../schemas/donors.schema.yaml")


###### Reference genome info ######
reference_name = config["reference"]["name"]
reference_fasta = config["reference"]["fasta"]

UT_subset = config["UT_subset"]

def get_mutect_and_delly_input(wildcards):
    """
    Mutect input. Read from sample sheet, return as a dict. 
    """
    input_dict = {"mutect": samples.loc[wildcards.pt]["mutect"],
                  "delly": samples.loc[wildcards.pt]["delly"]}
    return input_dict


def aggregate_preop_estimates(wildcards):
    """
    Pre-treatment/baseline cfDNA sample model estimates
    """
    preop_estimates = preop_units.preop_tf_path[preop_units.sample_ID == wildcards.pt].item()
    return preop_estimates