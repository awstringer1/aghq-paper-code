# Bayesian Inference using Adaptive Gauss-Hermite Quadrature

This repository contains the code used to implement the examples in *Higher Order Accurate Approximate Bayesian Inference using Adaptive Quadrature* by Blair Bilodeau, Alex Stringer, and Yanbo Tang (equal contribution authors), currently in preparation.

The examples all make use of the `aghq` package, available [on github](https://github.com/awstringer1/aghq) and to be available from `CRAN`.

## Code files

The following files may be run to recreate the examples from the paper:

- **Simulations (Section 5)**: file `11-simulations-final.R`.

- **Example 6.1**: file `12-disease-final.R`. The `TMB` function template `02_disease.cpp` must be compiled once beforehand, using the `compile()` call in the `.R` file.

- **Example 6.2**: file `13-astro-final.R`. The `TMB` function template `01_astro.cpp` must be compiled once beforehand, using the `compile()` call in the `.R` file.

- **Example 7.1**: file `10-zip-loaloa.R`. The `TMB` function template `10_loaloazip.cpp` must be compiled once beforehand, using the `compile()` call in the `.R` file.
