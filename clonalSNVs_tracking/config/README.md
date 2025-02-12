# General settings
To configure this workflow, modify ``config/config.yaml`` file according to your needs, which handles the input data configuration files and the reference genome details.

A complete list of configuration parameters with  human readable descriptions is available in the [config schema](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/workflow/schemas/config.schema.yaml).

### Reference genome
A name of the reference genome the data is aligned to and a path to the reference genome FASTA file needs to be provided in the config.

### Ensembl Variant Effect Predictor cache
For variant annotation, the workflow uses Ensembl Variant Effect Predictor (VEP)[1]. The transcript models file (cache) needs to be downloaded before using the workflow and the path to the cache needs to be specified in config.

See more how to download the cache from [VEP documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html).


# Input data

## Sample sheet
Patient data can be specified in the `config/samples.tsv` configuration file. For each patient, primary tumor and matched germline WGS data BAM files are required along with a variant set called with Mutect2[2]. The full schema for the `samples.tsv` file can be found [here](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/workflow/schemas/samples.schema.yaml).

Optionally, `samples.tsv` can also be used to define tumor purity and ploidy ranges for Sequenza, VAF cutoffs for candidate SNVs filtering and the cutoff for variants cluster assignments for MOBSTER. The default values for these parameters can be found in the [schema](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/workflow/schemas/samples.schema.yaml). 


## Units sheet
cfDNA data paths can be specified in the `config/units.tsv`. A number of cfDNA samples can be specified for each patient, by stating the corresponding patient `sample_ID` from `samples.tsv`, the unique `cfDNA_ID` of the sample, the cfDNA WGS BAM file path and the timepoint of the sample. Descriptions of the required parameters can be found in the [units schema](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/workflow/schemas/units.schema.yaml).



## References
[1] W. McLaren et al., “The Ensembl Variant Effect Predictor,” Genome Biol., vol. 17, no. 1, p. 122, Jun. 2016, doi: 10.1186/s13059-016-0974-4.
[2] K. Cibulskis et al., “Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples,” Nat. Biotechnol., vol. 31, no. 3, pp. 213–219, Mar. 2013, doi: 10.1038/nbt.2514.
