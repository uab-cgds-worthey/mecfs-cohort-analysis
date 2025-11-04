--------------------------------------------------------------------------------
CIBERSORTx Group Mode - OUTPUT FILES
--------------------------------------------------------------------------------

MAIN OUTPUT FILES
--------------------------------------------------------------------------------
CIBERSORTx Group Mode imputes representative cell type-specific expression 
profiles. In doing so, it generates a set of regression coefficients that 
represent the average expression value of each gene within each cell type across
the set (i.e., "group") of input mixture files.

The main result of CIBERSORTx Group Mode is a file for cell-type specific gene
expression profiles where genes have been filtered out using a threshold to
eliminate unreliably estimated genes for each cell type. This file is called
*_GEPs_Filtered.txt.

The "1" values in the expression matrix txt files are genes with insufficient
evidence of expression (these genes are either not expressed or have inadequate
statistical power to be imputed).

The NA values are genes that have inadequate statistical power to be imputed.

*_GEPs.txt is the file for cell-type specific gene expression profiles where no
filtering was done.

ADDITIONAL OUTPUT FILES
--------------------------------------------------------------------------------
- *_Fractions.txt: file enumerating the fractions of the different cell types in
bulks samples. 

The different statistics used for the filtering are saved in the files listed
below. Refer to Supplementary Note in Newman et al. (submitted) for further
details:
- *_GEPs_StdErrs.txt: analytically derived standard errors of the regression
coefficients.
- *GEPs_Pvals.txt: p-values used to determine the significance of the
regression coefficients.
- *GEPs_Qvals.txt: adjusted p-values (q-values) after multiple hypothesis
testing using the Benjamini-Hochberg method.
- *_GEPs_CV.txt: To further reduce confounding noise, genes were filtered based
on their geometric coefficient of variation (geometric c.v.), which are the
values listed in this file, calculated using the natural logarithm of
subsampled regression coefficients. The geometric c.v. were used to
determined the adaptive cell-type specific noise threshold.
- *_GEPs_ThresholdPlots.pdf: plot illustrating the adaptive noise threshold 
used for filtering.
- CIBERSORTxGEP_Weights.txt: the fractions of the different cell types after
merging them into major classes according to the merged classes file. 

If ground truth was given as input:
- *_CrossCorrelationMatrix.pdf: plots showing the correlation between estimated
genes and ground truth for each cell types. for all genes (GEP), and for the
signature matrix genes (SM). The corresponding *_CrossCorrelationMatrix.txt file
is also given as output.
- *_ScatterPlots.pdf: scatterplots showing the estimated expression values (non-
zero only) versus the observed expression values for the whole gene expression 
profile (GEP) and for the signature matrix genes (SM) after noise filtering.
- *_SM_GEPS_HeatMap.png: Heatmap illustrating the CiBERSORTx imputed gene
expression values for the signature matrix genes (y axis), compared to ground
truth.
- *GEPs_Stats.txt: set of benchmark statistics used to compare with ground truth.