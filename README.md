# MultinomialMutations

## Introduction

The R package `MultinomialMutations` was written to estimate a multinomial logistic regression model on the mutations (C to T, C to G, C to A, etc.) observed in cancer genomes, with genomic and epigenomic variables that describe the local environment, like nucleotide context, DNase1 hypersensitivity, etc, as regressors. Of course, other R packages for multinomial regression are available, but the huge data (hundreds or thousands of cancer genomes, with each genome consisting of approximately 3 billion positions) makes it necessary to implement a faster algorithm that takes advantage of the sparsity of the design matrix.

In the main function of the package, `fast_multinom`, the multinomial regression model with M categories is divided into M-1 binomial regression models. In each of them, parameters are estimated with the function `MatrixModels` from the contributed R-package [`MatrixModels`](https://cran.r-project.org/web/packages/MatrixModels/index.html). Next, the covariance matrix is computed for the joint model as described in Begg & Gray (1984). This is the core part of the R-package. It is a flexible function that is not restricted to the modeling of mutations, but it can be used for multinomial regression in a very general setting. The usage of the function mimics other regression functions like `lm` and `glm` and it also comes with a `predict` function.

For the analysis of the mutation data obtained from Fredriksson et al. (2014), data preparation and imputation functions have been added (although most of the data preparation is handled by -- *Guo's repo* --; `sumstats_orig`, `count_table_prep_multinom`, `mimp`, `add_counts`), as well as an example dataset (`cancermutations`). Furthermore, the deviance loss can be computed for variable selection by cross validation (`deviance_loss`) and estimates can be pooled after multiple imputation (`pooling`). Nested contrast matrices can be created with the functions `nested_sum_contrasts` and `nested_treatment_contrasts`. 

More details about the data and model are found in Bertl, J., Guo, Q. et al. (2017). The scripts and wrapper functions used for the analysis in the paper are included in the package in the folder `inst/Bertl_et_al_2017` and are explained in more detail below. The data will be made available. -- *figshare* --

### Installation

Download (clone) the R package into a folder, say `downloadfolder`, and then type in R:
```r
install.packages("downloadfolder/MultinomialMutations", repos=NULL)
```

Alternatively, install directly from github with the R package [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html):
```r
install.packages("devtools")
library(devtools)
install_github("MultinomialMutations/MultinomialMutations")
```
The package depends on the packages `psych, methods, Matrix, MatrixModels, data.table, reshape2, rmarkdown`. Note that the data preparation functions depend on `data.table` and `reshape2`, whereas the estimation function `fast_multinom` depends on `Matrix` and `MatrixModels`.


## Small usage example

For each genomic position, 4 outcomes are possible: no mutation (NO), transition (VA), transversion (VG) to a G:C basepair and transversion to a A:T basepair (VT). When the model includes the variable `strong` (dummy variable for G:C basepairs), this translates to a strand-symmetric mutation model. 

In this example, we estimate parameters to quantify the impact of the APOBEC mutation signature on the different mutation types. A cancer type specific intercept is included to take into account that the mutation rate varies among cancer types.

```r
# Load example data:
data(cancermutations)

# # The APOBEC mutation signature is only relevant for transitions (VA) and transversions to a G:C basepair (VG). Construct the corresponding subset of parameters for the 3 binomial models VA vs NO, VG vs NO and VT vs NO:
subs = matrix(T, ncol=3, nrow=4)
subs[3,2] = F

# Fit a multinomial logistic regression model to the mutation rate with the explanatory variables `strong` (C or G position), `apobec` (TpCpA or Tp**C**pT position) and a cancer type specific intercept.
fit = fast_multinom(cbind(NO, I, VA, VG) ~ strong + apobec + cancer_type, data = cancermutations, refLevel=1, loglik=T, predictions=T, VC=T, subsetmatrix=subs)

# Predictions on the data that was used for fitting (only available, because predictions=T in the function fast_multinom):
head(predict.fast_multinom(fit))
```



## Scripts used for the analysis in Bertl, J., Guo, Q. et al (2017)

-- *make data availabe on figshare and fill in here* --

## References 

### The package is used in: 

Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 

Juul, M.; Bertl, J.; Guo, Q.; Nielsen, M. M.; Świtnicki, M.; Hornshøj, H.; Madsen, T.; Hobolth, A. & Pedersen, J. S. Non-coding cancer driver candidates identified with a sample- and position-specific model of the somatic mutation rate. eLife, 2017, 6, e21778. https://doi.org/10.7554/eLife.21778

### Literature:

Begg, C. B. & Gray, R. Calculation of polychotomous logistic regression parameters using individualized regressions. Biometrika, 1984, 71, 11-18

Fredriksson, N. J.; Ny, L.; Nilsson, J. A. & Larsson, E. Systematic Analysis of noncoding somatic mutations and gene expression alterations across 14 tumor types. Nature Genetics, 2014, 46, 1258-1263