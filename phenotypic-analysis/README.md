# ME-CFS-Phenotypic-Comparison

A collection of diseases found in the Ramsey award ME/CFS cohort, phenotypes of ME/CFS and those diseases, and analysis
scripts for the comparison and visualization of phenotypes b/w diseases.

## Installation

### Requirements

- Anaconda3 or Mamba
- Git v2.0+

### Setup

Installation starts with fetching the Git repo and cloning it

```sh
git clone https://github.com/uab-cgds-worthey/mecfs-cohort-analysis.git
cd mecfs-cohort-analysis/phenotypic-analysis/
```

then the conda environment can built from the cloned repo's root directory like

```sh
# build environment using conda
conda env create -f config/conda-env.yml
# equivalent build command using mamba (which is generally much faster to build)
mamba env create -f config/conda-env.yml
```

the environemnt can be activated via

```sh
conda activate mecfs-pheno
```

## Phenotype Analyses

Phenotype analyses were completed using a combination of Python, R, and Jypter notebook. These analyses
leverage the relationship structure of the human phenotype ontology (HPO) for several quantitative
analyses described below.

### Avatar Disease Enrichment Ranking
<!-- markdown-link-check-disable -->
MECFS diagnosis criteria outlined in [Carruthers et. al (2011)](https://doi.org/10.1111/j.1365-2796.2011.02428.x)
<!-- markdown-link-check-enable -->
were [converted to HPO terms](data/mecfs-phenotypes.tsv) per diagnostic group and then used
to compute all possible minimal combinations that met Carruthers et. al (2011) diagnostic criteria, refered to as
MECFS phenotypic avatars or just avatars. These avatars could be used for similarity measurements with participants
in the Ramsey award cohort and also used for finding the most enriched (phenotypically overrepresented) diseases
across all avatars.

Deeper explination of the analytical process and computation used can be found in the
[avatar_enrichment Jupyter notebook](notebooks/avatar_enrichment.ipynb).

The resulting table of ranked enriched diseases from MECFS avatars will be saved to
[the results directory in an Excel file](data/results/Supp.%20Table%20-%20MECFS%20avatars%20and%20ranked%20disease%20enrichment%20table.xlsx).

### Cohort Disease Enrichment Ranking

Participant phenotypes were derived from symptoms noted by each participant in an in-take survey and stored per
participant in [a CSV file in this analysis](data/mecfs-samples/patient_hpo_summaries.csv). Like the avatars these
were used for finding the most enriched (phenotypically overrepresented) diseases across the cohort.

Deeper explination of the analytical process and computation used can be found in the
[avatar_enrichment Jupyter notebook](notebooks/avatar_enrichment.ipynb).

The resulting table of ranked enriched diseases from the cohort will be saved to
[the results directory in an Excel file](data/results/Supp.%20Table%20-%20MECFS%20participants%20ranked%20disease%20enrichment%20table.xlsx)

### Cohort Phenotypes

This analysis focuses on looking at participant phenotypic similarity unsing a binary encoding of phenotype terms.
