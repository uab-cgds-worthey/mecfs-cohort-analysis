# mecfs-dge-analysis

This repository hosts the differential gene expression analysis of a cohort of 29 patients who have been diagnosed with ME/CFS.

## Background

Myalgic Encephalomyelitis/Chronic Fatigue Syndrome (ME/CFS) is a chronic and debilitating illness affecting millions of Americans, characterized by severe fatigue, pain, flu-like symptoms, and cognitive issues. The cause of ME/CFS is not well understood, but there is evidence suggesting a genetic predisposition and dysregulation of the immune system leading to an overactive immune response.

## Our Hypothesis

We hypothesize that ME/CFS is a genetically inherited disease that results in disrupted metabolic regulation of the immune system,thus resulting in an inappropriately activated immune system.

## Rerunning the analysis

This project has been created using the workflowr template.

Analysis Rmarkdown files are located in the `analysis` folder.

To rerun the entire workflow, run the below in your R console':

```r
install.packages('workflowr')
library(workflowr)
wflow_build("analysis/*.Rmd")
wflow_view()
```

This will allow you to view the analysis website.

## Results

TBD.

## Authors

Shaurita D. Hutchins [:email:](mailto:sdhutchins@uab.edu) | Graduate Student
