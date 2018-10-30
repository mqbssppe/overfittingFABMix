# fabMix
 
## R code for *Overfitting Bayesian Mixtures of Factor Analyzers with an Unknown Number of Components* [arXiv:1701.04605 [stat.ME]](https://arxiv.org/abs/1701.04605).

Beta-version of the package is up:

1. [Download tarball](https://github.com/mqbssppe/overfittingFABMix/blob/master/fabMixPackage/version_2.0/fabMix_2.0.tar.gz)
2. [See Basic Usage](https://github.com/mqbssppe/overfittingFABMix/wiki/Basic-usage)
2. [Example with missing values](https://github.com/mqbssppe/overfittingFABMix/wiki/Example-with-missing-values)


__NOTE: all fabMix_versions < 4.2 are compatible only with Linux distributions.__

## Updates

### [September 2017] Version_2.0 added with support for **missing values** :waxing_gibbous_moon:

### [March 2018] Paper accepted to [Computational Statistics & Data Analysis](https://doi.org/10.1016/j.csda.2018.03.007) (please use version 2.0 of the software for the relevant scripts)

### [September 2018] Versions 3 and 4 added, which contain a plethora of new models and functionalities. 

Version 3.0 introduces new parsimonious models. 

Since version 4.0, the package is integrated with C++ code. 

Version 4.1 improves the Lambda output: all MCMC values are now exported to a single file (instead of multiple ones, as done in previous versions). 

### [October 2, 2018] fabMix version 4.2 added with Windows compatibility.

### [October 8, 2018] fabMix version 4.3 added (currently under development):

The `fabMix` function now features a new argument (`parallelModels`), allowing the user to run different models in parallel. This is combined with the pre-existing option to run heated chains in parallel, thus, parallelization is now implemented in __both model-level and chain-level__. 

### [October 30, 2018] fabMix 4.3 available on CRAN

[CRAN package](https://cran.r-project.org/package=fabMix)
