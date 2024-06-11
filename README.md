# mecfs-dge-analysis

This repository hosts the differential gene expression analysis of a cohort of 29 patients who have been diagnosed with ME/CFS.

## Background

Myalgic Encephalomyelitis/Chronic Fatigue Syndrome (ME/CFS) is a chronic and debilitating illness affecting millions of Americans, characterized by severe fatigue, pain, flu-like symptoms, and cognitive issues. The cause of ME/CFS is not well understood, but there is evidence suggesting a genetic predisposition and dysregulation of the immune system leading to an overactive immune response.

## Rerunning the analysis

This project has been created using the `workflowr` template.

The analysis and associated code is modified from [Dr.Gurpreet Kaur's](https://github.com/gurpreet-bioinfo) example: https://github.com/uab-cgds-worthey/bulk_rna-seq_dge

Analysis Rmarkdown files are located in the `analysis` folder.

To rerun the entire workflow, clone this repository and run the below in your R console. This project is using both `renv` and `workflowr` to improve reproducibility and allow results to be viewed in a web browser.

### Cloning this repository

```bash
git clone https://github.com/uab-cgds-worthey/mecfs-dge-analysis.git
cd mecfs-dge-analysis/
```

### Restoring the Environment

Open this project in RStudio. Upon opening, renv should automatically install.

If `renv` does not install automatically, run ```install.packages('renv')``` and 
proceed with the instructions below.

```r
# After cloning the project, open it in RStudio
# renv should automatically install
# restore the r environment
renv::restore()

# Rebuild the workflow R site to view all figures in the browser.
library(workflowr)
wflow_build("analysis/*.Rmd")
wflow_view()
```

## Results

Please review detailed results in our publication.

## Authors

Shaurita D. Hutchins [:email:](mailto:sdhutchins@uab.edu) | Graduate Student
