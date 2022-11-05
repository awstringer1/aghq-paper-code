# Bayesian Inference using Adaptive Gauss-Hermite Quadrature

This repository contains the code used to implement the examples in **Stochastic Convergence Rates and Applications of Adaptive Quadrature in Bayesian Inference** by Blair Bilodeau, Alex Stringer, and Yanbo Tang. A [preprint](https://arxiv.org/abs/2102.06801) and open access [journal version](https://www.tandfonline.com/doi/full/10.1080/01621459.2022.2141635) are available.

The examples all make use of the `aghq` package, available [on github](https://github.com/awstringer1/aghq) as well as `CRAN` (`install.packages("aghq")`). A [vignette](https://arxiv.org/abs/2101.04468) for that package is also available. The code currently here uses version `0.4.0` of `aghq`.

## Code files

The following files may be run to recreate the examples from the paper:

- **Section 4.1**: file `12-disease-final.R`.

- **Section 4.2**: file `13-astro-final.R`.

- **Section 5.2, S.7.1**: file `05-loaloazip.R`. 
  - Takes a few hours at the resolution shown in the paper; the default in the script
  produces lower-resolution spatial interpolations and runs in 5 -- 10 minutes or so.
  
  - The `MCMC` run took 66 hours on my server. By default the `domcmc` flag is set
  to `FALSE`.
  
- **Section S.7.2**: file `01-loaloazip-fixedintercept.R`.

- **Section S.6**: file `11-simulations-final.R`.

- **Section S.7.3**: file `06-loaloazipsims.R`
  - Takes a few hours on my server, minutes on my M1 MacBook Pro.

All files save results in a created subdirectory of `tempdir()`. To find them,
run `tempdir()` and go to that directory. The code will also print the location
once it's done running.

You will need to install other packages from `CRAN` as necessary within each script.
The only non-`CRAN` package is `ipoptr` which requires a working installation of
`IPOPT`, see [here](https://coin-or.github.io/Ipopt/INSTALL.html). This is laborious
and is only required for the astro example (Section 4.2). Everything else should
be pretty straightforward to run.
