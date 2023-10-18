![Bionano logo](images/Bionano-Logo.png?raw=true)


## Bionano - Clinical Affairs Virtual Karyotype Workflow with Snakemake
---
This repository contains a Snakemake workflow for the generation and visualization of virtual karyotypes.


## Overview

The workflow performs the following steps:

1. Run a virtual karyotype.
2. Generate plotting parameter files.
3. Visualize the virtual karyotype.
4. Rotate the ideogram images.
5. Merge rotated ideogram images.
6. Copy the input data to a separate location.
7. Package the results for easy sharing.

## Requirements

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* Conda (for environment management)
The specific Python, R, and tool dependencies are handled using Conda environments. For each rule, a conda environment is specified that will be automatically created and activated by Snakemake.

**Conda installation**

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh
#Follow-installation prompts
bash Miniconda3-py311_23.5.2-0-Linux-x86_64.sh
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


## Repository setup
---
1. Clone this repository and initialize OMKar submodule (VK github repository)
```
git clone https://joeyestabrook@bitbucket.org/bionanoclinicalaffairs/virtual_karyotype_snakemake.git
git submodule init
git submodule update
```
2. Modify the `config.yaml`: Specify the paths for the required resources (`centro`,`cytoband`,`resolved_cytoband`)

## Parameters
`centro`
* Description: Path to the file containing information on centromeres for the human genome (version hg38).

`cytoband`
* Description: Path to the cytoband file. This file provides cytogenetic banding pattern information, which is essential for generating the virtual karyotype plots.

`resolved_cytoband`
* Description: Path to the updated cytoband file, which contains 800 band resolution and is formatted to include gieStain and base cytoband ids.

`basedir`
* Description: Base directory for the virtual karyotype analyses. This might be the directory where all the analysis outputs will be saved or where some primary resources are located.

### Generating Input Files for `run_vk`` Rule
The input files required for the `run_vk` rule are derived from a de novo assembly. Here's how they are generated:

### De Novo Assembly:
Before using this workflow, a de novo assembly should be carried out. The output from this assembly will provide essential files that will be used as inputs for the `run_vk` rule.

### Directory Structure:
The relevant files from the de novo assembly can be found in the following directory structure:

* Smap File:
`output/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_filter_inversions.smap`

* Copy Number exp File:
`output/contigs/alignmolvref/copynumber/cnv_calls_exp.txt`

* Copy Number rcmap File:
`output/contigs/alignmolvref/copynumber/cnv_rcmap_exp.txt`

* Annotation Xmap File:
`output/contigs/annotation/exp_refineFinal1_merged.xmap`

Please ensure that these files are in the directory structure outlined below. Once these files are in place, you can run the `run_vk`.

3. Prepare your input samples: The workflow assumes input samples are placed under the samples directory with the structure:

```
samples/
    sample1/
        exp_refineFinal1_merged_filter_inversions.smap
        cnv_calls_exp.txt
        cnv_rcmap_exp.txt
        exp_refineFinal1_merged.xmap
    sample2/
    ...
```

## WORKFLOW OVERVIEW

![VK Dag](images/VK_dag.png?raw=true)


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
snakemake -p --verbose --use-conda -c 1 &> snakemake_workflow.log
```
Detach screen
```
ctrl-a ctrl-d
```

## Output
The primary output of this workflow includes:

* Processed virtual karyotype results in the vk_results directory for each sample.
* Packaged results for easy sharing in the packaged_results directory.

## Logging
Logs for each rule are saved in the `logs/scripts` directory. They can be used for troubleshooting and monitoring the progress of the workflo