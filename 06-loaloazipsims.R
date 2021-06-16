### Simulated data from the ZIP loaloa binomial model ###
# MCMC was/is diverging but aghq is working. To assess this in more detail,
# we simulate from the model and compare results to the truth

set.seed(75489)

## Load the data ##
library(tidyverse)
library(TMB)
precompile()
library(geostatsp)
data(loaloa,package = "geostatsp")
library(aghq)
library(tmbstan)

savestamp <- "20210505-v1"
globalpath <- normalizePath(tempdir(),winslash='/')
plotpath <- normalizePath(file.path(globalpath,"loaloazip"),winslash='/',mustWork = FALSE)
if (!dir.exists(plotpath)) dir.create(plotpath)
savepath <- plotpath

file.copy(
  normalizePath(system.file('extsrc/05_loaloazip.cpp',package='aghq'),winslash='/'),
  globalpath
)

# Compile TMB template-- only need to do once
compile(normalizePath(file.path(globalpath,"05_loaloazip.cpp"),winslash='/'))
dyn.load(normalizePath(dynlib(file.path(globalpath,"05_loaloazip")),winslash='/'))

## Setup ----

## Prepare the "inner" model ##

# Design matrices
Amat <- Diagonal(nrow(loaloa))

Xmat <- cbind(rep(1,nrow(Amat)))
# Design matrix: zip model and risk model are the same
design <- bdiag(
  # ZIP
  cbind(
    Amat,
    Xmat
  ),
  # Risk
  cbind(
    Amat,
    Xmat
  )
)

# Response
y <- loaloa@data$y
N <- loaloa@data$N

## Dimensions
n <- nrow(Xmat) # Number of obs
p <- ncol(Xmat) * 2 # Number of betas
m <- ncol(Amat) * 2 # Number of spatial points
Wd <- ncol(design) # Number of total params
# Check
stopifnot(Wd == m + p)

## Prior distributions ##
# Use the same prior for both sets of Matern params
priorparams <- list(
  sigma_u = 1,
  sigma_alpha = .025,
  rho_u = 2e05,
  rho_alpha = .975,
  beta_prec = 1 / 2^2, # SD of 2
  nu = 1
)


# PC Prior for kappa,tau
maternconstants <- list()
maternconstants$d <- 2 # Dimension of spatial field, fixed
maternconstants$nu <- 1 # Shape param, fixed
get_kappa <- function(sigma,rho) sqrt(8*maternconstants$nu)/rho
get_tau <- function(sigma,rho) sigma * get_kappa(sigma,rho)^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d/2) * (4*pi)^(maternconstants$d/2) / gamma(maternconstants$nu))
get_sigma <- function(kappa,tau) tau / (kappa^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d/2) * (4*pi)^(maternconstants$d/2) / gamma(maternconstants$nu)))
get_rho <- function(kappa,tau) sqrt(8*maternconstants$nu) / kappa


simulate_spatial_fields <- function(U,
                                    theta,
                                    pointsdata,
                                    resolution = list(nrow = 100,ncol = 100)) {
  # U: matrix of samples, each column is a sample
  # theta: data.frame of theta values
  # Draw from U*|U
  
  # Compute matrix of var, range, shape
  modpar <- cbind(
    var = get_sigma(exp(theta$theta1),exp(theta$theta2))^2,
    range = get_rho(exp(theta$theta1),exp(theta$theta2)),
    shape = maternconstants$nu
  )
  
  fielddat <- pointsdata
  fielddat@data <- as.data.frame(U)
  
  geostatsp::RFsimulate(
    model = modpar,
    data = fielddat,
    x = raster(fielddat,nrows = resolution$nrow,ncols = resolution$ncol)
  )
}

# Fit the original model
sigma_u <- 1
sigma_alpha <- .025
rho_u <- 2e05
rho_alpha <- .975
beta_prec <- .001

# PC Prior for kappa,tau
maternconstants <- list()
maternconstants$d <- 2 # Dimension of spatial field, fixed
maternconstants$nu <- 1 # Shape param, fixed
get_kappa <- function(sigma,rho) sqrt(8*maternconstants$nu)/rho
get_tau <- function(sigma,rho) sigma * get_kappa(sigma,rho)^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d/2) * (4*pi)^(maternconstants$d/2) / gamma(maternconstants$nu))
get_sigma <- function(kappa,tau) tau / (kappa^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d/2) * (4*pi)^(maternconstants$d/2) / gamma(maternconstants$nu)))
get_rho <- function(kappa,tau) sqrt(8*maternconstants$nu) / kappa

datlist <- list(
  y = y,
  N = N,
  design = design,
  nu = maternconstants$nu,
  rho_u = rho_u,
  rho_alpha = rho_alpha,
  sigma_u = sigma_u,
  sigma_alpha = sigma_alpha,
  D = raster::pointDistance(loaloa,lonlat = FALSE),
  betaprec = beta_prec
)
startingsig <- 1
startingrho <- 4.22*1e04
paraminit <- list(
  W = rnorm(ncol(design)),
  logkappa = log(get_kappa(startingsig,startingrho)),
  logtau = log(get_tau(startingsig,startingrho))
)

ff <- MakeADFun(data = datlist,
                parameters = paraminit,
                random = "W",
                DLL = "05_loaloazip",
                ADreport = FALSE,
                silent = TRUE)

loaloazipquad <- aghq::marginal_laplace_tmb(
  ff,
  3,
  startingvalue = c(paraminit$logkappa,paraminit$logtau)
)



## Simulate Data ----
ilogit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log( x / (1-x) )
sample_zip_binomial <- function(phi,p,Ni) {
  # phi: probability zero is haphazard
  out <- rbinom(length(Ni),Ni,p)
  out[runif(length(phi)) > phi] <- 0
  out
}

origsamps <- sample_marginal(loaloazipquad,1e03)
origWmean <- apply(origsamps$samps,1,mean)
origeta <- design %*% origWmean
origzipprobs <- ilogit(origeta[1:190])
origrisk <- ilogit(origeta[191:380])


simulate_model <- function(pars) {
  # randomly sample ntosim points with replacement
  # and then simulate zero-inflated binomial response
  # ntosim <- 190
  itr <- pars['itr']
  ntosim <- pars['ntosim']
  n <- nrow(loaloa)
  # Get a vector of how many times to sample each village
  # https://stackoverflow.com/questions/24845909/generate-n-random-integers-that-sum-to-m-in-r
  xx <- sort(sample(1:(ntosim - 1*(n-1)),n-1,replace = TRUE),decreasing = TRUE)
  xx <- c(xx[1:(n-2)] - xx[2:(n-1)] + 1,xx[n-1])
  xx <- c(xx,ntosim - sum(xx))
  stopifnot(length(xx) == n)
  stopifnot(sum(xx) == ntosim)
  stopifnot(all(xx >= 1))
  
  
  # idx <- sort(sample(1:n,ntosim,replace = TRUE))
  idx <- rep(1:n,times = xx)
  idxu <- unique(idx)
  
  # idx <- 1:190
  rowidx <- c(idx,n+idx)
  colidx <- c(idxu,n+1,n+1+idxu,2*(n+1))
  
  # pull the correct probs
  sample_design <- design[rowidx,colidx]
  sampleeta <- as.numeric(sample_design %*% origWmean[colidx])
  samplezipprobs <- ilogit(sampleeta[1:ntosim])
  samplerisk <- ilogit(sampleeta[(ntosim+1):(2*ntosim)])
  
  sampledata <- loaloa[idx, ]
  sampledatau <- loaloa[idxu, ]
  
  
  # simulate the data
  simdataframe <- SpatialPointsDataFrame(sampledata,
                                         data = data.frame(
                                           N = sampledata$N,
                                           y = sample_zip_binomial(samplezipprobs,samplerisk,sampledata$N)
                                         ))
  
  # jitter the points to remove overlap
  # simdataframe@coords <- simdataframe@coords + matrix(rnorm(2*nrow(simdataframe),sd = 1000),ncol = 2)
  
  # Fit the model
  datlist <- list(
    y = simdataframe$y,
    N = simdataframe$N,
    design = sample_design,
    nu = maternconstants$nu,
    rho_u = priorparams$rho_u,
    rho_alpha = priorparams$rho_alpha,
    sigma_u = priorparams$sigma_u,
    sigma_alpha = priorparams$sigma_alpha,
    D = raster::pointDistance(sampledatau,lonlat = FALSE),
    betaprec = priorparams$beta_prec
  )
  
  # set.seed(456443)
  

  
  loaloazipquad_sim <- list()
  class(loaloazipquad_sim) <- "try-error"
  inneritr <- 1
  while("try-error" %in% class(loaloazipquad_sim) & inneritr < 5) {
    cat("ntosim =",ntosim,", itr =",itr,", inneritr =",inneritr,".\n")
    paraminit <- list(
      W = rnorm(ncol(sample_design),sd = .1),
      logkappa = log(get_kappa(startingsig,startingrho)),
      logtau = log(get_tau(startingsig,startingrho))
    )
    
    ff <- MakeADFun(data = datlist,
                    parameters = paraminit,
                    random = "W",
                    DLL = "05_loaloazip",
                    ADreport = FALSE,
                    silent = TRUE)
    tm <- Sys.time()
    loaloazipquad_sim <- try(aghq::marginal_laplace_tmb(
      ff,
      3,
      startingvalue = c(paraminit$logkappa,paraminit$logtau)
    ),silent = TRUE)
    modtime <- difftime(Sys.time(),tm,units = 'secs')
    inneritr <- inneritr + 1
  }
  failed <- "try-error" %in% class(loaloazipquad_sim)
  
  filepathout <- file.path(plotpath,paste0("rmse-sim-table",savestamp,".csv"))
  if (failed) {
    rmse <- -1
    covr <- -1
    rmse_suit <- -1
    covr_suit <- -1
    rmse_inc <- -1
    covr_inc <- -1
    out <- data.frame(n = ntosim,iter = itr,rmse=rmse,covr = covr,rmse_suit=rmse_suit,covr_suit = covr_suit,rmse_inc=rmse_inc,covr_inc = covr_inc,time = modtime)
    readr::write_csv(out,filepathout,append = file.exists(filepathout))
    return(out)
  }

  loazippostsamples <- sample_marginal(loaloazipquad_sim,1e04)
  
  estW <- apply(loazippostsamples$samps,1,mean)
  Wlower <- apply(loazippostsamples$samps,1,quantile,probs = .025)
  Wupper <- apply(loazippostsamples$samps,1,quantile,probs = .975)
  
  # Note: use original design matrix, not sample design matrix. Works
  # because of ordering of the sample
  EE <- design %*% loazippostsamples$samps
  # Rows 1:190 are suitability; 191:380 (conditional) incidence
  probs <- apply(EE,1,ilogit)

  suitprobsmean <- apply(probs[ ,1:190],2,mean)
  suitprobslower <- apply(probs[ ,1:190],2,quantile,probs = .025)
  suitprobsupper <- apply(probs[ ,1:190],2,quantile,probs = .975)
  
  incprobsmean <- apply(probs[ ,191:380],2,mean)
  incprobslower <- apply(probs[ ,191:380],2,quantile,probs = .025)
  incprobsupper <- apply(probs[ ,191:380],2,quantile,probs = .975)
  
  
  # RMSE
  rmse <- sqrt(mean( (estW - origWmean[colidx])^2 ))
  rmse_suit <- sqrt(mean( (suitprobsmean - origzipprobs)^2 ))
  rmse_inc <- sqrt(mean( (incprobsmean - origrisk)^2 ))
  
  # Coverage
  covr <- mean(origWmean[colidx] <= Wupper & origWmean[colidx] >= Wlower)
  covr_suit <- mean(origzipprobs <= suitprobsupper & origzipprobs >= suitprobslower)
  covr_inc <- mean(origrisk <= incprobsupper & origrisk >= incprobslower)
  
  out <- data.frame(n = ntosim,iter = itr,rmse=rmse,covr = covr,rmse_suit=rmse_suit,covr_suit = covr_suit,rmse_inc=rmse_inc,covr_inc = covr_inc,time = modtime)
  readr::write_csv(out,filepathout,append = file.exists(filepathout))
  out
}

simstodo <- expand.grid(1:50,c(200,500,1000,5000,10000))

simlist <- list()
for (i in 1:nrow(simstodo)) simlist[[i]] <- c("itr" = simstodo[i,1],"ntosim" = simstodo[i,2])

library(parallel)
options(mc.cores = parallel::detectCores())
sims <- mclapply(simlist,simulate_model)
save(sims,file = file.path(savepath,paste0("simresults",savestamp,".RData")))

## Make the boxplots ----
simresults <- readr::read_csv(file.path(plotpath,paste0("rmse-sim-table",savestamp,".csv")))

# rmse = -1 indicates failed run; how many?
simresults %>% filter(rmse == -1) %>% group_by(n) %>% summarize(howmany = n())

# 2021/04/25 seems ok
# n howmany
# <dbl>   <int>
# 1   200       6
# 2   500       1

DEFAULTTHEME <- theme_classic() +
  theme(text = element_text(size = 28))

simresultstoplot <- simresults %>%
  mutate(nn = factor(n,levels = sort(unique(simresults$n)))) %>%
  filter(rmse != -1)

# RSMELABS <- labs(x = "Simulated Sample Size",y = "RMSE")
# COVRLABS <- labs(x = "Simulated Sample Size",y = "Coverage")
RSMELABS <- labs(x = "",y = "")
COVRLABS <- labs(x = "",y = "")

# W
rmseWboxplot <- simresultstoplot %>%
  ggplot(aes(x = nn,y = rmse)) +
  DEFAULTTHEME +
  geom_boxplot() +
  RSMELABS

covWboxplot <- simresultstoplot %>%
  ggplot(aes(x = nn,y = covr)) +
  scale_y_continuous(breaks = seq(.9,1,by=.01)) +
  coord_cartesian(ylim = c(.9,1)) +
  DEFAULTTHEME +
  geom_boxplot() +
  geom_hline(yintercept = .95,linetype = "dashed",colour = "gray") +
  COVRLABS

# Suitability
rmsesuitboxplot <- simresultstoplot %>%
  ggplot(aes(x = nn,y = rmse_suit)) +
  DEFAULTTHEME +
  geom_boxplot() +
  RSMELABS

covsuitboxplot <- simresultstoplot %>%
  ggplot(aes(x = nn,y = covr_suit)) +
  scale_y_continuous(breaks = seq(.9,1,by=.01)) +
  coord_cartesian(ylim = c(.9,1)) +
  DEFAULTTHEME +
  geom_boxplot() +
  geom_hline(yintercept = .95,linetype = "dashed",colour = "gray") +
  COVRLABS

# Incidence
rmseincboxplot <- simresultstoplot %>%
  ggplot(aes(x = nn,y = rmse_inc)) +
  DEFAULTTHEME +
  geom_boxplot() +
  RSMELABS

covincboxplot <- simresultstoplot %>%
  ggplot(aes(x = nn,y = covr_inc)) +
  scale_y_continuous(breaks = seq(.9,1,by=.01)) +
  coord_cartesian(ylim = c(.9,1)) +
  DEFAULTTHEME +
  geom_boxplot() +
  geom_hline(yintercept = .95,linetype = "dashed",colour = "gray") +
  COVRLABS

# Save them
ggsave(plot = rmseWboxplot,file = file.path(plotpath,paste0("loaloa-W-rmse-boxplot-",savestamp,".pdf")),width = 7,height = 7)
ggsave(plot = covWboxplot,file = file.path(plotpath,paste0("loaloa-W-cov-boxplot-",savestamp,".pdf")),width = 7,height = 7)
ggsave(plot = rmsesuitboxplot,file = file.path(plotpath,paste0("loaloa-suit-rmse-boxplot-",savestamp,".pdf")),width = 7,height = 7)
ggsave(plot = covsuitboxplot,file = file.path(plotpath,paste0("loaloa-suit-cov-boxplot-",savestamp,".pdf")),width = 7,height = 7)
ggsave(plot = rmseincboxplot,file = file.path(plotpath,paste0("loaloa-inc-rmse-boxplot-",savestamp,".pdf")),width = 7,height = 7)
ggsave(plot = covincboxplot,file = file.path(plotpath,paste0("loaloa-inc-cov-boxplot-",savestamp,".pdf")),width = 7,height = 7)

