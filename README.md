![Bionano logo](images/Bionano-Logo.png?raw=true)


## Bionano - Clinical Affairs Virtual Karyotype Snakemake Workflow
---
This pipeline is designed to visualize a virtual karyotype using Denovo assembly output files.

## Repository setup
---
```
git clone https://joeyestabrook@bitbucket.org/bionanoclinicalaffairs/virtual_karyotype_snakemake.git
git submodule init
git submodule update
```

## SETUP
---
Current implementation of this workflow utilizes the following packages
```
conda (4.12.0)
mamba (0.15.3)
Python (3.10.6)
```


**Conda installation**

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
#Follow-installation prompts
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh
```

**Mamba installation**

```
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
This will install snakemake into an isolated software environment, that has to be activated with

```
$ conda activate snakemake
$ snakemake --help
```


## RUNNING WORKFLOW
---

Do a dry-run of snakemake to ensure proper execution prior to launching the workflow.

```
screen
```

```
$(snakemake) snakemake -np --verbose
```
Execute snakemake workflow
```
snakemake -p --verbose --use-conda -c 1
```
Detach screen
```
ctrl-a ctrl-d
```