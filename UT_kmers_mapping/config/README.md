# General settings
To configure this workflow, modify ``config/config.yaml`` file according to your needs, which handles the input data configuration files and the reference genome details.

A complete list of configuration parameters with  human readable descriptions is available in the [config schema](https://github.com/BesenbacherLab/ctDNAmer/tree/main/UT_kmers_mapping/workflow/schemas/config.schema.yaml).

### Reference genome
A name of the reference genome the data is aligned to and a path to the reference genome FASTA file needs to be provided in the config.


# Input data

## Sample sheet
Patient data can be specified in the `config/samples.tsv` configuration file. For each patient the UT set in text file format created with the ctDNAmer workloq is needed. The full schema for the `samples.tsv` file can be found [here](https://github.com/BesenbacherLab/ctDNAmer/tree/main/UT_kmers_mapping/workflow/schemas/samples.schema.yaml).

