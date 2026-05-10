# Nextstrain build for ANDV

## About 

Based on CCHF build - will update the readme more later.

This repository provides scripts for the phylogenetic analysis of the CCHF virus. CCHF is a multi-segmented RNA-virus with three segments (L, M and S). Phylogenetic analysis is performed on each segment using the references; [NC_005301.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_005301.3) for segment L, [NC_005300.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_005300.2) for segment M and [NC_005302.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_005302.1) for segment S. This repository creates nextstrain builds with or without linking the segments.

When linking is specified trees are built using samples (isolates) which have all three segments. These samples are given the display name `isolate name/collection date/country` in each tree, allowing them to be visualized using a pairwise tanglegram.

In the case where linking is not specified segment trees are built individually and we annotate lineages (sometimes called "genotypes") of the Segment S using the annotation defined in 

__Serena A. Carroll, Brian H. Bird, Pierre E. Rollin, Stuart T. Nichol,
Ancient common ancestry of Crimean-Congo hemorrhagic fever virus__ (2010), [link](https://doi.org/10.1016/j.ympev.2010.01.006).

The resulting, non-linked trees can be used in a nextclade dataset.


## Installation

To run the workflow using TreeKnit you need to set up an environment containing [Nextstrain](https://docs.nextstrain.org/en/latest/install.html) and Julia with [TreeKnit](https://github.com/PierreBarrat/TreeKnit.jl) and TreeTools. The workflow can then be build from within this environment.

### WARNING
TreeKnit currently only runs on julia 1.7 - this does not exist micromamba for arm, 
```
micromamba create -n julia17_trial -c conda-forge julia=1.7 --platform osx-64
```
creates an environment with julia 1.7, however this isn't compatible with the other dependencies. I installed julia 1.7.3 via dmg and then compiled the treeknit package to run these scripts.

### Local Development

For local development activate the micromamba environment (this includes all required packages except treeknit) using:

```
    micromamba create -f environment.yml --platform=linux/amd64
    micromamba activate cchfv
```

You can format your snakemake files using

```
snakefmt ${file}
```

You can visualize the workflow using:

```
snakemake -s ${file} --dag | dot -Tpng > dag.png
```

## Run the Workflow

There are three different snakefile versions provided within this repository.

The standard snakefile called `Snakefile` contains the whole workflow including a resolving step using TreeKnit. To use TreeKnit it is necessary to call a Julia script and thus Julia has to be installed in the environment.

If julia is not available or TreeKnit should be excluded from the workflow the snakefile `Snakefile_without_treeknit` can be used alternatively, by using the file option of snakemake:

```
-s FILE
```

The third snakefile `Snakefile_without_linking` creates the phylogenetic trees of all segments separately without ensuring that the sequences for each segment are chosen from the same set of samples, making it harder to create tanglegrams of the trees.
