# Bayesian Inference using Adaptive Gauss-Hermite Quadrature

This repository contains the code used to implement the examples in *Higher Order Accurate Approximate Bayesian Inference using Adaptive Quadrature* by Blair Bilodeau, Alex Stringer, and Yanbo Tang (equal contribution authors), currently in preparation.

The examples all make use of the `aghq` package, available [on github](https://github.com/awstringer1/aghq) and to be available from `CRAN`.

## Code files

The following files may be run to recreate the examples from the paper:

- **Section 4.1**: file `12-disease-final.R`.

- **Section 4.2**: file `13-astro-final.R`.

- **Section 5.2**: file `05-loaloazip.R`. 
  - Takes a few hours at the resolution shown in the paper; the default in the script
  produces lower-resolution spatial interpolations and runs in 5 -- 10 minutes or so.
  
  - The `MCMC` run took 66 hours on my server. By default the `domcmc` flag is set
  to `FALSE`.

- **Section S.5**: file `11-simulations-final.R`.

- **Section S.6**: file `06-loaloazipsims.R`
  - Takes a few hours.

All files save results in a created subdirectory of `tempdir()`. To find them,
run `tempdir()` and go to that directory.

You will need to install other packages from `CRAN` as necessary within each script.
The only non-`CRAN` package is `ipoptr` which requires a working installation of
`IPOPT`, see [here](https://coin-or.github.io/Ipopt/INSTALL.html). This is laborious
and is only required for the astro example (Section 4.2). Everything else should
be pretty straightforward to run.