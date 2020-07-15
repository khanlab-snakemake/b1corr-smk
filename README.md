# Snakemake workflow: b1corr-smk
Snakemake workflow for analyses as performed in ['Effects of MP2RAGE B1+ sensitivity on inter-site T1 reproducibility and morphometry at 7T'](https://doi.org/10.1101/2020.02.13.947382)

## Authors

* Roy AM Haast @royhaast 

## Usage

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yml` to configure the workflow execution, and `participants.tsv` to specify your subjects.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda --use-singularity -n

Execute the workflow locally via

    snakemake --use-conda --use-singularity --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --use-singularity --cluster qsub --jobs 100

or

    snakemake --use-conda --use-singularity --drmaa --jobs 100


If you are using Compute Canada, you can use the [cc-slurm](https://github.com/khanlab/cc-slurm) profile, which submits jobs and takes care of requesting the correct resources per job (including GPUs). Once it is set-up with cookiecutter, run:

    snakemake --profile cc-slurm

Or, with [neuroglia-helpers](https://github.com/khanlab/neuroglia-helpers) can get a 1-GPU, 8-core, 32gb node and run locally there. First, get a node with a GPU (default 8-core, 32gb, 3 hour limit):

    regularInteractive -g
    
Then, run:

    snakemake --use-conda --use-singularity --cores 8 --resources gpu=1 mem=32000


See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

