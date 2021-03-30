# Bayesian Inference using Adaptive Gauss-Hermite Quadrature

This repository contains the code used to implement the examples in **Stochastic Convergence Rates and Applications of Adaptive Quadrature in Bayesian Inference** by Blair Bilodeau, Alex Stringer, and Yanbo Tang (equal contribution authors). The paper has been submitted to a journal and a [preprint](https://arxiv.org/abs/2102.06801) is available.

The examples all make use of the `aghq` package, available [on github](https://github.com/awstringer1/aghq) as well as `CRAN` (`install.packages("aghq")`). A [vignette](https://arxiv.org/abs/2101.04468) for that package is also available.

**Important note**: the results in the [preprint](https://arxiv.org/abs/2102.06801) all reproduce exactly when using the `CRAN` version (0.1.0) of the package. The package is under active development and the development version may not always give the same results. If you want to reproduce the results as seen in the [preprint](https://arxiv.org/abs/2102.06801), use the `CRAN` version.

## Code files

The following files may be run to recreate the examples from the paper:

- **Simulations (Section 5)**: file `11-simulations-final.R`.

- **Example 6.1**: file `12-disease-final.R`. The `TMB` function template `02_disease.cpp` must be compiled once beforehand, using the `compile()` call in the `.R` file. The saved `MCMC` results can be downloaded [here](https://drive.google.com/file/d/1_ogYlDHsmQYTviCLqdwiCrxOhLmHoTJ7/view?usp=sharing).

- **Example 6.2**: file `13-astro-final.R`. The `TMB` function template `01_astro.cpp` must be compiled once beforehand, using the `compile()` call in the `.R` file.

- **Example 7.1**: file `10-zip-loaloa.R`. The `TMB` function template `10_loaloazip.cpp` must be compiled once beforehand, using the `compile()` call in the `.R` file.
