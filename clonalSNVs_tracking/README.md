# cfDNA circulating tumor fraction estimation based on clonal SNVs allele frequencies

This workflow estimates circulating tumor fraction from algined cfDNA WGS data based on clonal SNVs allele frequencies. The clonal SNVs are identified from the aligned primary tumor WGS data, quality filtered and then tracked in the cfDNA samples. 

## Description
The method is built as a customizable [snakemake workflow](https://f1000research.com/articles/10-33)[1]. 

This workflow first identifies clonal SNVs from primary tumor data and then tracks these SNVs in the cfDNA by following these main steps:
1) Input Mutect2[2] calls filtering with vcftools to retain only SNV variants (indels removed) with PASS FILTER flag, located on autosomes or allosomes. 
2) Variant annotation with Ensembl[Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)[3]. Variants in coding regions with a "high" impact annotation are marked as (potential) drivers. 
3) Calling Copy Number Alterations with the [Sequenza toolset](https://sequenzatools.bitbucket.io/#/home)[4].
4) Quality control of the SNVs, CNAs, tumor purity and ploidy estimates using [CNAqc](https://caravagnalab.github.io/CNAqc/)[5].
5) Selecting SNVs in diploid heterozygous (1:1) copy number states.
6) Removing outlier SNVs based on VAF distribution. 
7) Finding clonal mutations (seperation from noise and passenger variants at lower frequencies) with [MOBSTER](https://caravagnalab.github.io/mobster/)[6]. 
8) Quality filtering the clonal SNVs based on primary tumor alignment information (looking at sequencing and alignment quality at both read and position level).
9) Tracking the clonal SNVs in cfDNA (includes quality filters on the cfDNA alignment).
10) Calculating the mean allele frequency of clonal SNV in cfDNA.

Detailed descripton of the workflow can be found in the Methods section of [ctDNA-mers publication](https://www.biorxiv.org/).

## Requirements

### Software
The workflow requires Snakemake 8.0.0 or above and uses conda for package management. 

### Input data

Patient data can be specified in the ``config/samples.tsv`` configuration file. For each patient, primary tumor and matched germline WGS data BAM files are required along with a variant set called with Mutect2[2]. 
``Samples.tsv`` can also be used to define tumor purity and ploidy ranges for Sequenza and VAF cutoffs for candidate SNVs filtering. The default values for these parameters can be found in the [samples schema](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/workflow/schemas/samples.schema.yaml). 

The ``config/config.yaml`` file specifies the name and FASTA file path of the reference genome where the samples are aligned to. 
For variant annotation, the workflow uses Ensembl Variant Effect Predictor (VEP)[3]. The transcript models file (cache) needs to be downloaded before using the workflow and the path to the cache needs to be specified in ``config.yaml``.

## Usage

### Step 1: Install snakemake

Snakemake is best to be installed via the [Mamba package manager](https://github.com/mamba-org/mamba) (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). For other options see [here](https://github.com/mamba-org/mamba).

Given that Mamba is installed, run:

```
mamba create -c conda-forge -c bioconda --name snakemake 'snakemake>=8'
```

to install Snakemake in an isolated environment. If you need to use conda instead of mamba, `--conda-frontend conda` flag needs to be added to the snakemake commands given below. 

Activate the environment via: 

```
conda activate snakemake
```

### Step 2: Clone this repo

Download and extract the parent repository: 

```
git clone https://github.com/carmenoroperv/ctDNA_mers.git && cd ctDNA_mers/clonalSNVs_tracking
```

### Step 3: Configure workflow 

#### Workflow confifuration
To specify the parameters for running the workdlow and the sample paths, modify the configuration files [`config.yaml`](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/config/config.yaml), [`samples.tsv`](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/config/samples.tsv) and [`units.tsv`](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/config/units.tsv) according to your needs, following the explanations provided [here](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/config).

### Step 4: Run workflow 

#### Cluster exection

For cluster execution of the workflow, the [snakemake slurm executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) needs to be installed with `pip install snakemake-executor-plugin-slurm`. If the slurm plugin is not installed, the `-e` flag needs to be specified for the snakemake commands listed below. 

The specifics for cluster execution should be defined in the workflow profile configuration file. An example workflow profile for slurm is provided [here](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/workflow/profiles/default/config.yaml). To use the example profile, adjust the snakemake command line parameters to your needs. Importantly, a cluster account is specified in the example profile as an environment variable. To set the account name as an environment variable run `export ACCOUNT_NAME=<your_account_name>` or modify the profile config file to include your account name directly.


After you have activated the conda environment with snakemake, installed the slurm executor plugin and set the account name as an environment variable, you can test the workflow remote execution by performing a dry-run:

```
snakemake -n
```

To run the workflow for a new data set, use the `--directory` flag that specifies the path to the directory where the pipeline will be executed. The target directory needs to include a config folder with `config.yaml` and `samples.tsv` files, which specify the ctDNA-mers parameters and paths to sample files that will be used during execution. You can execute the workflow with: 

```
snakemake --directory "path/to/new/directory/"
```

The workflow profile that specifies the details for the cluster execution will be still automatically detected from the pipeline directory ([`workflow/profiles/default/config.yaml`](https://github.com/carmenoroperv/ctDNA_mers/tree/main/clonalSNVs_tracking/workflow/profiles/default/config.yaml)) even when the execution directory is changed. If you want to specify a new cluster execution profile as well, use the `--workflow-profile` flag: 

```
snakemake --workflow-profile "path/to/workflow_profile/config.yaml"
```

For further options for local, cluster and cloud execution, see the snakemake [docs](https://snakemake.readthedocs.io/).


## References
[1] F. Mölder et al., “Sustainable data analysis with Snakemake,” F1000Research, vol. 10, p. 33, Apr. 2021, doi: 10.12688/f1000research.29032.2.
[2] K. Cibulskis et al., “Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples,” Nat. Biotechnol., vol. 31, no. 3, pp. 213–219, Mar. 2013, doi: 10.1038/nbt.2514.
[3] W. McLaren et al., “The Ensembl Variant Effect Predictor,” Genome Biol., vol. 17, no. 1, p. 122, Jun. 2016, doi: 10.1186/s13059-016-0974-4.
[4] F. Favero et al., “Sequenza: allele-specific copy number and mutation profiles from tumor sequencing data,” Ann. Oncol., vol. 26, no. 1, pp. 64–70, Jan. 2015, doi: 10.1093/annonc/mdu479.
[5] A. Antonello et al., “Computational validation of clonal and subclonal copy number alterations from bulk tumor sequencing using CNAqc,” Genome Biol., vol. 25, no. 1, p. 38, Jan. 2024, doi: 10.1186/s13059-024-03170-5.
[6] G. Caravagna et al., “Subclonal reconstruction of tumors by using machine learning and population genetics,” Nat. Genet., vol. 52, no. 9, pp. 898–907, Sep. 2020, doi: 10.1038/s41588-020-0675-5.