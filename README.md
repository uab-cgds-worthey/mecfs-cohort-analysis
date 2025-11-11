# ME/CFS N-of-1

<!-- markdown-link-check-disable -->

[![R lintr](https://github.com/uab-cgds-worthey/mecfs-cohort-analysis/actions/workflows/lintr.yaml/badge.svg)](https://github.com/uab-cgds-worthey/mecfs-cohort-analysis/actions/workflows/lintr.yaml)
[![DGE Results Deployment](https://github.com/uab-cgds-worthey/mecfs-cohort-analysis/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/uab-cgds-worthey/mecfs-cohort-analysis/actions/workflows/pages/pages-build-deployment)
[![Linting - Markdown,Shell](https://github.com/uab-cgds-worthey/mecfs-cohort-analysis/actions/workflows/linting.yml/badge.svg)](https://github.com/uab-cgds-worthey/mecfs-cohort-analysis/actions/workflows/linting.yml)

<!-- markdown-link-check-enable -->

This repository contains multiple complementary analyses of the Ramsey award ME/CFS cohort, including differential
gene expression analysis and phenotypic comparison studies.

## Background

Myalgic Encephalomyelitis/Chronic Fatigue Syndrome (ME/CFS) is a chronic and debilitating illness affecting millions of
individuals worldwide. It is characterized by severe fatigue, pain, flu-like symptoms, and cognitive issues. The cause
of ME/CFS is not well understood, but evidence suggests a genetic predisposition and dysregulation of the immune
system leading to an overactive immune response.

This repository hosts analyses that explore:

- Differential gene expression patterns in ME/CFS patients
- Phenotypic overlap between ME/CFS and other diseases
- Patient phenotypic similarity and clustering
- Disease enrichment analyses

## Repository Structure

- `dge-analysis/` - Differential gene expression analysis using R
- `phenotypic-analysis/` - Phenotypic comparison analysis
- `variant-findings-plot/` - Variant findings plotting

## Installation

### Requirements

**For Differential Gene Expression Analysis:**

- R ([v4.5.1](https://cran.r-project.org/bin/macosx/big-sur-x86_64/base/R-4.5.1-x86_64.pkg) or later)
- RStudio (v2025.05.1+513 or later)

**For Phenotypic Comparison Analysis and Variant Findings plot:**

- Anaconda3 or Mamba

**Common Requirements:**

- Git v2.0+

### Setup

Installation starts with fetching the Git repo and cloning it:

```bash
git clone https://github.com/uab-cgds-worthey/mecfs-cohort-analysis.git
cd mecfs-cohort-analysis/
```

For detailed setup instructions for each analysis, see the README files in their respective directories:

- **Differential Gene Expression Analysis:** See [`dge-analysis/README.md`](dge-analysis/README.md) for detailed
  R/renv setup instructions
- **Phenotypic Comparison Analysis:** See the [phenotypic comparison documentation](phenotypic-analysis/README.md) for
  conda/mamba environment setup instructions
- **Variant Findings plot:** See the [variant findings documentation](variant-findings-plot/README.md) for conda/mamba
  environment setup instructions

## Analyses

### Differential Gene Expression Analysis

The differential gene expression analysis examines transcriptional changes in a cohort of 23 patients diagnosed with
ME/CFS using bulk RNA-seq. The analysis identifies genes that are differentially expressed between
affected and unaffected conditions, and explores subgroup-specific patterns.

Located in [`dge-analysis/`](dge-analysis/).

For detailed documentation including setup, input data, workflow instructions, and results, see [`dge-analysis/README.md`](dge-analysis/README.md).

#### Results: Differential Gene Expression

- Results and figures can be found in the `dge-analysis/docs/` folder
- Viewable online at: [https://uab-cgds-worthey.github.io/mecfs-cohort-analysis/](https://uab-cgds-worthey.github.io/mecfs-cohort-analysis/)

### Phenotypic Comparison Analysis

A collection of diseases found in the Ramsey award ME/CFS cohort, phenotypes of ME/CFS and those diseases, and
analysis scripts for the comparison and visualization of phenotypes between diseases.

Located in [`phenotypic-analysis/`](phenotypic-analysis/)

For detailed documentation including setup, tools, notebooks, and analysis instructions, see
[`phenotypic-analysis/README.md`](phenotypic-analysis/README.md).

#### Results: Phenotypic Comparison

- Results and figures can be found in [`phenotypic-analysis/data/results/`](phenotypic-analysis/data/results/)
- Carruthers ME/CFS diagnositic creteria as HPO terms is defined in
  [phenotypic-analysis/data/mecfs-phenotypes.tsv](phenotypic-analysis/data/mecfs-phenotypes.tsv)

### Variant findings plot

Organized metadata and primary variant findings organized in a variant format expected by the pyoncoprint
Python library for generating variant oncoplots. The figure generation process is described in the
[plotting Jupyter notebook](variant-findings-plot/notebooks/mecfs_oncoprint.ipynb)

## Authors

- Shaurita D. Hutchins [:email:](mailto:sdhutchins@uab.edu) | PhD Candidate
- Brandon M. Wilk [:email:](mailto:bwilk777@uab.edu) | PhD Candidate

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).
