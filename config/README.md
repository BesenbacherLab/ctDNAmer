# General settings
To configure the workflow, modify ``config/config.yaml``, which handles the input data configuration files and ctDNA-mers parameter settings.

A complete list of configuration parameters with human readable descriptions is available in the [config schema](https://github.com/carmenoroperv/ctDNA_mers/tree/main/workflow/schemas/config.schema.yaml).

# Input data configuration files

## Sample sheet
Patient data can be specified in the `config/samples.tsv`. For each patient, primary tumor and matched germline WGS data are required. The full schema for the `samples.tsv` file can be found [here](https://github.com/carmenoroperv/ctDNA_mers/tree/main/workflow/schemas/samples.schema.yaml). The input data can be provided in FASTQ or BAM format. 

## Units sheet
cfDNA data paths can be specified in the `config/units.tsv`. A number of cfDNA samples can be specified for each patient, by stating the corresponding patient `sample_ID` from `samples.tsv`, the unique `cfDNA_ID` of the sample, the cfDNA WGS data files and the timepoint of the cfDNA sample. The input data can be provided in FASTQ or BAM format. Descriptions of the required parameters can be found in the [units schema](https://github.com/carmenoroperv/ctDNA_mers/tree/main/workflow/schemas/units.schema.yaml).

## Germline union samples sheet
Samples that will be used to create the germline union set can be specified in `config/samples_glu.tsv`. These germline samples will be combined to a union k-mer set and applied during germline subtraction for each patient in the target patient cohort.
For each germline sample a unique sample_ID, type of sample (e.g. buffycoat or cfDNA) and the path to the data files need to be specified. The data can be provided either in FASTQ or BAM format. Descriptions of the required parameters can be found in the [schema](https://github.com/carmenoroperv/ctDNA_mers/tree/main/workflow/schemas/samples_glu.schema.yaml).

## Donors sheet
Unmatched cfDNA samples for empirical noise distribution estimation can be provided in `config/donors.tsv`. For each sample, a unique donor_ID and cfDNA data file paths need to be provided. 
The data can be provided either in FASTQ or BAM format. Descriptions of the required parameters can be found in the [schema](https://github.com/carmenoroperv/ctDNA_mers/tree/main/workflow/schemas/donors.schema.yaml).


# Count filter test settings
To configure the setup for the count-filter test, modify ``config/config_count_filter_test.yaml``, which handles the input data configuration files and ctDNA-mers parameter settings. The input fields of the count-filter test configuration are identical to the ctDNA-mers tool. The patients and cfDNA samples can be specified separately for the count-filter test in `config/samples_count_filter_test.tsv` and `config/units_count_filter_test.tsv`. The same germline union and donor samples for empirical noise estimation are applied for the count-filter test as for the ctDNA-mers tool. This default behaviour can be adjusted by creating additional sample configuration files and specifying the paths to these in ``config/config_count_filter_test.yaml``. A complete list of configuration parameters with  human readable descriptions is available in the [config schema](https://github.com/carmenoroperv/ctDNA_mers/tree/main/workflow/schemas/config_count_filter_test.schema.yaml).

## Input data

### Sample sheet
Patient data can be specified in the `config/samples_count_filter_test.tsv`. The full schema for the `samples_count_filter_test.tsv` file can be found [here](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/workflow/schemas/samples_count_filter_test.schema.yaml).

### Units sheet
cfDNA data paths can be specified in the `config/units_count_filter_test.tsv`. A cfDNA sample can be specified by stating the corresponding patient `sample_ID` from `samples_count_filter_test.tsv`, the unique `cfDNA_ID` of the sample, the cfDNA WGS data files and the timepoint of the sample. Descriptions of the required parameters can be found in the [schema](https://github.com/carmenoroperv/ctDNA_mers/tree/main/workflow/schemas/units_count_filter_test.schema.yaml).

At least two cfDNA samples are required per patient for testing: a ctDNA-positive pre-treatment/baseline cfDNA sample and a ctDNA-negative/post-treatment cfDNA sample. The count filter test runs TF estimation for unique tumor sets of different sizes and the minimum required unique tumor set size can be determined based on the ctDNA-positive and -negative samples TF estimates difference. 


## References
[1] W. McLaren et al., “The Ensembl Variant Effect Predictor,” Genome Biol., vol. 17, no. 1, p. 122, Jun. 2016, doi: 10.1186/s13059-016-0974-4.
[2] K. Cibulskis et al., “Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples,” Nat. Biotechnol., vol. 31, no. 3, pp. 213–219, Mar. 2013, doi: 10.1038/nbt.2514.
