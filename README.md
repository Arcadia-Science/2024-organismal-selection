# 2024-organismal-selection

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

Note: Analysis repo names should be prefixed with the year (e.g., `2024-noveltree-analysis`). This prefix can be changed at time of publication if appropriate.

## Purpose

This repository contains code for performing analyses related to the pub "Leveraging evolution to identify novel organismal models of human biology". Code for proteome curation, phylogenomic inference, and molecular conservation calculations can be found HERE [replace with link].

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/).

After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), you can now build the environment. Because the conservation analysis depends on several R-packages not distributed through conda, as well as several software that must be locally compiled from source, you must take two additional steps before building the environment. First, you must edit the [environment YAML file](./envs/aa_stats_mv_dists.yml), uncommenting the C/C++ compilers that are appropriate for your operating system. This section of the environment YAML file is shown below.

```
dependencies: # Uncomment below based on your OS
  #- gcc_linux-64 # Uncomment if running on Linux
  #- gxx_linux-64 # Uncomment if running on Linux, GCC C++ Compiler
  #- clang_osx-64 # Uncomment if running on macOS
  #- clangxx_osx-64 # Uncomment if running on macOS, Clang C++ Compiler
```

Second, you must run an additional [build script](./code/build_remaining_env_aa_stats_mv_dists.sh) after creating and activating the new conda environment with the appropriate compilers installed. Below, we provide code to carry out the whole process (after modifying the environment YAML file).

```{bash}
# Create the environment and activate it
mamba env create -n aa_stats_mv_dists --file envs/aa_stats_mv_dists.yml
conda activate aa_stats_mv_dists

# Install the remaining software within this conda environment:
bash ./scripts/build_remaining_env_aa_stats_mv_dists.sh

# RYAN: DESCRIBE YOUR CONDA ENV HERE.
```

## Data

Before proceeding with any (re)analysis, first download the NovelTree run outputs from Zenodo HERE and decompress the outputs # TODO

```
# !!!!!  DOWNLOAD FROM ZENODO  !!!!! #

# !!!!! DECOMPRESS THE TARBALL !!!!! #
```

This directory, `2024-organismal-selection-zenodo`, contains the following:

- `run_configurations/noveltree-model-euks-samplesheet.csv` - the samplesheet for our snakemake preprocessing workflow to filter and preprocess species proteomes prior to analysis with NovelTree.
- `run_configurations/euk_preprocess_samplesheet.tsv` & `run_configurations/noveltree-model-euks-parameterfile.json` - the NovelTree sample and parameter files used to run NovelTree.
- `preprocessed_proteomes.tar.gz` - a compressed tarball containing the preprocessed proteomes used by our NovelTree run.
- `results-noveltree-model-euks.tar.gz` - a compressed tarball containing all outputs generated by our NovelTree run.
- `aa-summary-stats.tar.gz` - a compressed tarball containing all AA summary statistics generated by `code/genefam_aa_summaries.py`.
- `gf-aa-multivar-distances.tar.gz` - a compressed tarball containing all result files produced by `code/calc_protein_mv_distances.R`.

TODO: RYAN, CONTINUE THIS WITH THE CONTENTS YOU ADD

## Usage

With the NovelTree run outputs downloaded and extracted into the base directory of this repository, we now proceed by calling the script `code/genefam_aa_summaries.py`. This bash script calculates for each protein sequence within each gene family, summaries of AA composition, as well as AA physical properties. All code below assumes that you have downloaded and extracted the directory `2024-organismal-selection-zenodo/` from this pubs correspoding Zenodo repository.

```{bash}
# Ensure we are calling this script within the correct conda environment
conda activate aa_stats_mv_dists

# Set the MSA directory to variable
msa_dir="2024-organismal-selection-zenodo/results-noveltree-model-euks/witch_alignments/original_alignments/"

# Now, run the script to calculate the physicochemical properties of each protein using ProtParam
python code/genefam_aa_summaries.py -t 10 $msa_dir
```

This will create a new directory called "aa-summary-stats/" that contains the calculated AA properties for each protein, and summarized for each gene famly. With these protein properties curated, we can now proceed with the calculation of pairwise multivariate distances between proteins within each gene family.

```{bash}
Rscript code/calc_protein_mv_distances.R
```

Briefly, this script:

1. Reads in the species tree from the NovelTree run results and time-calibrates it using a species tree containing these species obtained from timetree.org
2. Reads in species metadata from the NovelTree samplesheet and copy number information
3. Reads in the gene family trees and protein properties calculated by `code/genefam_aa_summaries.py`, retaining only those gene families that contain human proteins, and then for each gene family, it:
   - Time-calibrates the gene family trees so branch lengths reflect time, rather than the extent of sequence divergence.
   - Uses this tree to transforms the AA physical properties such that we correct for phylogenetic non-independence between proteins
   - Calculate multivariate (mahalanobis) distances between proteins

TODO: RYAN, ADD IN YOUR CODE USAGE BELOW.

```{bash}
conda activate RYANS_ENVIRONMENT

MORE THINGS
```

## Overview

### Description of the folder structure

### Methods

TODO: Include a brief, step-wise overview of analyses performed.

> Example:
>
> 1. Download scripts using `download.sh`.
> 2. Preprocess using `./preprocessing.sh -a data/`
> 3. Run analysis script using `analysis.Rscript`
> 4. Generate figures using `pub/make_figures.R`.

### Compute Specifications

TODO: Describe what compute resources were used to develop and run the analysis. For example, you could list the operating system, number of cores, RAM, and storage space. You should log any major changes to the compute specifications here as they happen.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).

---

## For Developers

This section contains information for developers who are working off of this template. Please adjust or edit this section as appropriate when you're ready to share your repo.

### GitHub templates

This template uses GitHub templates to provide checklists when making new pull requests. These templates are stored in the [.github/](./.github/) directory.

### `.gitignore`

This template uses a `.gitignore` file to prevent certain files from being committed to the repository.

### Linting

This template automates linting using GitHub Actions and the [`lintr` linter](https://cran.r-project.org/web/packages/lintr/vignettes/lintr.html). When you push changes to your repository, GitHub will automatically run the linter and report any errors, blocking merges until they are resolved.
