# Unique tumor k-mers alignment to the reference venome

This workflow aligns unique tumor (UT) k-mers to the human reference genome. 

## Requirements

### Software
The workflow requires Snakemake 8.0.0 or above and uses conda for package management. 

### Input data

Patient data can be specified in the ``config/samples.tsv`` configuration file. 

The ``config/config.yaml`` file specifies the name and FASTA file path of the reference genome where the samples are aligned to. 

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
git clone https://github.com/BesenbacherLab/ctDNAmer.git && cd ctDNAmer/UT_kmers_mapping
```

### Step 3: Configure workflow 

#### Workflow confifuration
To specify the parameters for running the workdlow and the sample paths, modify the configuration files [`config.yaml`](https://github.com/BesenbacherLab/ctDNAmer/tree/main/UT_kmers_mapping/config/config.yaml) and [`samples.tsv`](https://github.com/BesenbacherLab/ctDNAmer/tree/main/UT_kmers_mapping/config/samples.tsv) according to your needs, following the explanations provided [here](https://github.com/BesenbacherLab/ctDNAmer/tree/main/UT_kmers_mapping/config).

### Step 4: Run workflow 

#### Cluster exection

For cluster execution of the workflow, the [snakemake slurm executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) needs to be installed with `pip install snakemake-executor-plugin-slurm`. If the slurm plugin is not installed, the `-e` flag needs to be specified for the snakemake commands listed below. 

The specifics for cluster execution should be defined in the workflow profile configuration file. An example workflow profile for slurm is provided [here](https://github.com/BesenbacherLab/ctDNAmer/tree/main/UT_kmers_mapping/workflow/profiles/default/config.yaml). To use the example profile, adjust the snakemake command line parameters to your needs. Importantly, a cluster account is specified in the example profile as an environment variable. To set the account name as an environment variable run `export ACCOUNT_NAME=<your_account_name>` or modify the profile config file to include your account name directly.


After you have activated the conda environment with snakemake, installed the slurm executor plugin and set the account name as an environment variable, you can test the workflow remote execution by performing a dry-run:

```
snakemake -n
```

To run the workflow for a new data set, use the `--directory` flag that specifies the path to the directory where the pipeline will be executed. The target directory needs to include a config folder with `config.yaml` and `samples.tsv` files, which specify the ctDNA-mers parameters and paths to sample files that will be used during execution. You can execute the workflow with: 

```
snakemake --directory "path/to/new/directory/"
```

The workflow profile that specifies the details for the cluster execution will be still automatically detected from the pipeline directory ([`workflow/profiles/default/config.yaml`](https://github.com/BesenbacherLab/ctDNAmer/tree/main/UT_kmers_mapping/workflow/profiles/default/config.yaml)) even when the execution directory is changed. If you want to specify a new cluster execution profile as well, use the `--workflow-profile` flag: 

```
snakemake --workflow-profile "path/to/workflow_profile/config.yaml"
```

For further options for local, cluster and cloud execution, see the snakemake [docs](https://snakemake.readthedocs.io/).

