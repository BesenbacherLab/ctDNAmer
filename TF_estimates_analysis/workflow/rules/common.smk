import pandas as pd
from pathlib import Path
import os
import sys
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("8.0.0")


###### Config file and sample sheets #####
configfile: "config/config_TF_estimates_analysis.yaml"
#validate(config, schema="../schemas/config.schema.yaml")

cohort = config["cohort"]
pref = config["pref"]

cfDNA_mean_plot = config["cfDNA_mean_plot"]