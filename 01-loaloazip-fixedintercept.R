### Loa Loa Zip Model ###
# Fix the intercepts for MCMC at the AGHQ estimates, to get MCMC to run.

## Load the data ##
library(TMB)
precompile()
library(geostatsp)
data(loaloa,package = "geostatsp")
library(aghq)
library(tmbstan)

library(parallel)
options(mc.cores = parallel::detectCores())

library(ggplot2)


# Set the resolution for the spatial interpolations.
# The results shown in the paper use:
reslist <- list(nrow = 200,ncol = 400)
# but this takes a couple hours. Here I set:
# reslist <- list(nrow = 50,ncol = 100)
# which should take only a few minutes
# Make sure you install the RandomFields package:
# install.packages('RandomFields')
# as this will speed up geostatsp::RFsimulate()
# Note that these considerations are specific to this example's
# post-processing and not related to the aghq method or package.

savestamp <- "20220116-v1"
globalpath <- normalizePath(tempdir(),winslash='/')
plotpath <- normalizePath(file.path(globalpath,"loaloazip"),winslash='/',mustWork = FALSE)
# globalpath <- "~/phd/projects/best-friends-gang/normalizing-constant"
# plotpath <- normalizePath(file.path(globalpath,"figures/revision"),winslash='/',mustWork = FALSE)
if (!dir.exists(plotpath)) dir.create(plotpath)
savepath <- plotpath
# savepath <- "~/data/loaloazip"
tmbpath <- tempdir()


file.copy(
  normalizePath(system.file('extsrc/05_loaloazip.cpp',package='aghq'),winslash='/'),
  tmbpath
)
file.copy(
  normalizePath(system.file('extsrc/05_loaloazip_fixedparam.cpp',package='aghq'),winslash='/'),
  tmbpath
)


# Compile TMB template-- only need to do once
compile(normalizePath(file.path(tmbpath,"05_loaloazip.cpp"),winslash='/'))
dyn.load(normalizePath(dynlib(file.path(tmbpath,"05_loaloazip")),winslash='/'))
compile(normalizePath(file.path(tmbpath,"05_loaloazip_fixedparam.cpp"),winslash='/'))
dyn.load(normalizePath(dynlib(file.path(tmbpath,"05_loaloazip_fixedparam")),winslash='/'))

# Flags, which analysis to do?
doaghq <- TRUE
domcmc <- TRUE
dopostsamplingaghq <- TRUE
dopostsamplingmcmc <- TRUE

# Initialize time variables
aghqtime <- 0
mcmctime <- 0
aghqsimtime <- 0
mcmcsimtime <- 0

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
sigma_u <- 1
sigma_alpha <- .025
rho_u <- 2e05
rho_alpha <- .975

# PC Prior for kappa,tau
maternconstants <- list()
maternconstants$d <- 2 # Dimension of spatial field, fixed
maternconstants$nu <- 1 # Shape param, fixed
get_kappa <- function(sigma,rho) sqrt(8*maternconstants$nu)/rho
get_tau <- function(sigma,rho) sigma * get_kappa(sigma,rho)^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d/2) * (4*pi)^(maternconstants$d/2) / gamma(maternconstants$nu))
get_sigma <- function(kappa,tau) tau / (kappa^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d/2) * (4*pi)^(maternconstants$d/2) / gamma(maternconstants$nu)))
get_rho <- function(kappa,tau) sqrt(8*maternconstants$nu) / kappa

# Precision for betas

beta_prec <- .001

## Log Posterior ----

startingsig <- 1
startingrho <- 4.22*1e04

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
# NOTE: for some initial values of W, TMB's inner optimization seems to fail
# This was tried over a bunch of random intializations and most worked, and all
# gave the same optimum. But this is why we set the seed here and use a random start.
set.seed(4564)
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

# New TMB template- create and test
betazifix <- 2.912111
betariskfix <- -1.984655

paraminitnew <- with(paraminit,list(
  Uzi = W[1:190],
  betazi = betazifix,
  Urisk = W[192:381],
  betarisk = betariskfix,
  logkappa = log(get_kappa(startingsig,startingrho)),
  logtau = log(get_tau(startingsig,startingrho))
))

datlistnew <- within(datlist,{
  A <- as(as(Amat,'dgCMatrix'),'dgTMatrix')
  X <- Xmat
})
datlistnew$design <- NULL
ffnew <- MakeADFun(data = datlistnew,
                parameters = paraminitnew,
                # random = c("Uzi","Urisk","betarisk",'betazi'),
                map = list(betazi = factor(NA),betarisk = factor(NA)),
                DLL = "05_loaloazip_fixedparam",
                ADreport = FALSE,
                silent = TRUE)

# stopifnot(as.numeric(with(ff,fn(par))) == as.numeric(with(ffnew,fn(par))))

if (doaghq) {
  tm <- Sys.time()
  cat("Doing AGHQ, time = ",format(tm),"\n")
  loaloazipquad <- aghq::marginal_laplace_tmb(
    ff,
    7,
    startingvalue = c(paraminit$logkappa,paraminit$logtau)
  )
  aghqtime <- difftime(Sys.time(),tm,units = 'secs')
  save(loaloazipquad,file = file.path(savepath,paste0("loaloazipquad",savestamp,".RData")))
  cat("AGHQ took: ",format(aghqtime),"\n")
}

## MCMC ----
# Do MCMC by FIXING the betas at their AGHQ estimates.
if (domcmc){
  paraminitmcmc <- paraminitnew
  paraminitmcmc$betarisk <- ff$env$last.par.best[382]
  paraminitmcmc$betazi <- ff$env$last.par.best[191]

  ffmcmc <- MakeADFun(data = datlistnew,
                     parameters = paraminitmcmc,
                     map = list(betazi = factor(NA),betarisk = factor(NA)),
                     DLL = "05_loaloazip_fixedparam",
                     ADreport = FALSE,
                     silent = TRUE)

  tm <- Sys.time()
  cat("Doing MCMC, time = ",format(tm),"\n")
  loaloazipmcmc <- tmbstan(
    ffmcmc,
    chains = 8,
    cores = 8,
    iter = 1e04,
    warmup = 1e03,
    init = Reduce(c,paraminit),
    seed = 124698,
    algorithm = "NUTS"
  )
  mcmctime <- difftime(Sys.time(),tm,units = 'secs')
  save(loaloazipmcmc,file = file.path(savepath,paste0("loaloazipmcmc",savestamp,".RData")))
  cat("MCMC took: ",format(mcmctime),"\n")
}

simulate_spatial_fields <- function(U,
                                    theta,
                                    pointsdata,
                                    resolution = list(nrow = 100,ncol = 100)) {
  # U: matrix of samples, each column is a sample
  # theta: data.frame of theta values
  # Draw from U*|U

  # Compute matrix of var, range, shape
  modpar <- cbind(
    var = get_sigma(exp(theta$logkappa),exp(theta$logtau))^2,
    range = get_rho(exp(theta$logkappa),exp(theta$logtau)),
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


## Post samples ----
if (FALSE) {
  # Load if in interactive mode (and you want to)
  load(file.path(savepath,paste0("loaloazipquad",savestamp,".RData")))
  load(file.path(savepath,paste0("loaloazipmcmc",savestamp,".RData")))
}

if (dopostsamplingaghq) {
  cat("Doing AGHQ field simulations, time = ",format(tm),"\n")
  loazippostsamples <- sample_marginal(loaloazipquad,100)
  # Extract out the U and V
  Uzi_idx <- which(rownames(loazippostsamples$samps)=="Uzi")
  betazi_idx <- which(rownames(loazippostsamples$samps)=="betazi")
  Urisk_idx <- which(rownames(loazippostsamples$samps)=="Urisk")
  betarisk_idx <- which(rownames(loazippostsamples$samps)=="betarisk")

  postU <- loazippostsamples$samps[Uzi_idx, ]
  postV <- loazippostsamples$samps[Urisk_idx, ]

  postBetazi <- loazippostsamples$samps[betazi_idx, ]
  postBetarisk <- loazippostsamples$samps[betarisk_idx, ]

  tm <- Sys.time()
  fieldbrickzip <- simulate_spatial_fields(
    postU,
    loazippostsamples$theta,
    loaloa,
    resolution = reslist
  )
  fieldbrickrisk <- simulate_spatial_fields(
    postV,
    loazippostsamples$theta,
    loaloa,
    resolution = reslist
  )

  fieldbrickzip_withint <- fieldbrickzip + postBetazi
  fieldbrickrisk_withint <- fieldbrickrisk + postBetarisk

  simfieldsmeanzip <- mean(1 / (1 + exp(-fieldbrickzip_withint)))

  simfieldsmeanrisk <- mean(
    (1 / (1 + exp(-fieldbrickzip_withint))) *
      (1 / (1 + exp(-fieldbrickrisk_withint)))
  )
  aghqsimtime <- difftime(Sys.time(),tm,units = 'secs')
  raster::writeRaster(fieldbrickzip,file.path(savepath,paste0("fieldbrickzipaghq",savestamp,".grd")),overwrite = TRUE)
  raster::writeRaster(fieldbrickrisk,file.path(savepath,paste0("fieldbrickzipaghq",savestamp,".grd")),overwrite = TRUE)

  cameroonBorderLL = raster::getData("GADM", country=c('CMR'), level=2)
  nigeriaBorderLL = raster::getData("GADM", country=c('NGA'), level=2)
  cameroonBorder = spTransform(cameroonBorderLL, projection(loaloa))
  nigeriaBorder = spTransform(nigeriaBorderLL, projection(loaloa))
  cameroonBorderouter <- rgeos::gUnaryUnion(cameroonBorder)
  nigeriaBorderouter <- rgeos::gUnaryUnion(nigeriaBorder)


  fullborder <- raster::bind(cameroonBorder,nigeriaBorder)
  fullborderouter <- raster::bind(cameroonBorderouter,nigeriaBorderouter)

  fullborder <- crop(fullborder,loaloa)
  fullborderouter <- crop(fullborderouter,loaloa)

  plot_loaloa <- function(plotraster,breaks) {
    predcols <- mapmisc::colourScale(
      plotraster,
      breaks = breaks,
      style = "fixed",
      col = "Spectral",
      rev = TRUE
    )

    plotraster <- mask(plotraster,fullborderouter)

    mapmisc::map.new(loaloa)
    plot(plotraster,
         col = predcols$col,
         breaks = predcols$breaks,
         legend=FALSE, add=TRUE)
    points(loaloa,pch = 4)
    plot(fullborder,add = TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
    plot(fullborderouter,add = TRUE)
    mapmisc::legendBreaks('right', predcols, cex=1, bty='o',inset = .05)
  }

  br <- c(0,.05,.1,.15,.2,.25,.3,.4,.5,.6,1)
  brzero <- c(.35,.5,.6,.7,.8,.9,.91,.92,.93,.94,.95,1)

  pdf(file.path(plotpath,paste0("loaloa-zip-postmean.pdf")),width=7,height=7)
  plot_loaloa(simfieldsmeanzip,brzero)
  dev.off()
  pdf(file.path(plotpath,paste0("loaloa-risk-postmean.pdf")),width=7,height=7)
  plot_loaloa(simfieldsmeanrisk,br)
  dev.off()
  cat("AGHQ simulations took: ",format(aghqsimtime),"\n")
}

## Post samples ----
if (dopostsamplingmcmc) {
  cat("Doing MCMC field simulations, time = ",format(tm),"\n")

  samps <- as.data.frame(loaloazipmcmc)
  Uzi_idx <- grep('Uzi',colnames(samps))
  Urisk_idx <- grep("Urisk",colnames(samps))
  theta_idx <- grep("log",colnames(samps))

  # Only do 100
  idx <- sample.int(nrow(samps),100)

  postU <- t(samps[idx,Urisk_idx])
  postV <- t(samps[idx,Uzi_idx])

  # Fixed betas
  postBetazi <- paraminitmcmc$betazi
  postBetarisk <- paraminitmcmc$betarisk

  theta <- samps[idx,theta_idx]

  tm <- Sys.time()
  fieldbrickzip <- simulate_spatial_fields(
    postU,
    theta,
    loaloa,
    resolution = reslist
  )
  fieldbrickrisk <- simulate_spatial_fields(
    postV,
    theta,
    loaloa,
    resolution = reslist
  )
  cat(difftime(Sys.time(),tm,units='secs'),'seconds')


  fieldbrickzip_withint <- fieldbrickzip + postBetazi
  fieldbrickrisk_withint <- fieldbrickrisk + postBetarisk

  simfieldsmeanzip <- mean(1 / (1 + exp(-fieldbrickzip_withint)))

  simfieldsmeanrisk <- mean(
    (1 / (1 + exp(-fieldbrickzip_withint))) *
      (1 / (1 + exp(-fieldbrickrisk_withint)))
  )
  mcmcsimtime <- difftime(Sys.time(),tm,units = 'secs')
  raster::writeRaster(fieldbrickzip,file.path(savepath,paste0("fieldbrickzipmcmc",savestamp,".grd")),overwrite = TRUE)
  raster::writeRaster(fieldbrickrisk,file.path(savepath,paste0("fieldbrickzipmcmc",savestamp,".grd")),overwrite = TRUE)

  cameroonBorderLL = raster::getData("GADM", country=c('CMR'), level=2)
  nigeriaBorderLL = raster::getData("GADM", country=c('NGA'), level=2)
  cameroonBorder = spTransform(cameroonBorderLL, projection(loaloa))
  nigeriaBorder = spTransform(nigeriaBorderLL, projection(loaloa))
  cameroonBorderouter <- rgeos::gUnaryUnion(cameroonBorder)
  nigeriaBorderouter <- rgeos::gUnaryUnion(nigeriaBorder)


  fullborder <- raster::bind(cameroonBorder,nigeriaBorder)
  fullborderouter <- raster::bind(cameroonBorderouter,nigeriaBorderouter)

  fullborder <- crop(fullborder,loaloa)
  fullborderouter <- crop(fullborderouter,loaloa)

  plot_loaloa <- function(plotraster,breaks) {
    predcols <- mapmisc::colourScale(
      plotraster,
      breaks = breaks,
      style = "fixed",
      col = "Spectral",
      rev = TRUE
    )

    plotraster <- mask(plotraster,fullborderouter)

    mapmisc::map.new(loaloa)
    plot(plotraster,
         col = predcols$col,
         breaks = predcols$breaks,
         legend=FALSE, add=TRUE)
    points(loaloa,pch = 4)
    plot(fullborder,add = TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
    plot(fullborderouter,add = TRUE)
    mapmisc::legendBreaks('right', predcols, cex=1, bty='o',inset = .05)
  }

  br <- c(0,.05,.1,.15,.2,.25,.3,.4,.5,.6,1)
  brzero <- c(.35,.5,.6,.7,.8,.9,.91,.92,.93,.94,.95,1)

  pdf(file.path(plotpath,paste0("loaloa-zip-postmean-mcmc.pdf")),width=7,height=7)
  plot_loaloa(simfieldsmeanzip,brzero)
  dev.off()
  pdf(file.path(plotpath,paste0("loaloa-risk-postmean-mcmc.pdf")),width=7,height=7)
  plot_loaloa(simfieldsmeanrisk,br)
  dev.off()
  cat("MCMC simulations took: ",format(mcmcsimtime),"\n")
}

# Write the timing table
timingtable <- data.frame(
  task = c("AGHQ","MCMC","AGHQSim","MCMCSim"),
  time = c(aghqtime,mcmctime,aghqsimtime,mcmcsimtime)
)
readr::write_csv(timingtable,file.path(savepath,paste0("timing-table-fixedparam",savestamp,".csv")))


## Get the comparison between MCMC and AGHQ posteriors at the sampled sites ##
nsamp <- nrow(samps)
aghqpostsamps <- sample_marginal(loaloazipquad,nsamp)
aghqpostsamplist <- split(aghqpostsamps$samps,1:nrow(aghqpostsamps$samps))
aghqpostsamplist <- aghqpostsamplist[grep("U",rownames(aghqpostsamps$samps))]
mcmcsamplist <- as.list(samps)
mcmcsamplist <- mcmcsamplist[grep("U",names(mcmcsamplist))]

get_ks <- function(x,y) as.numeric(ks.test(x,y)$statistic)
ks_stats <- mcmapply(get_ks,aghqpostsamplist,mcmcsamplist)
ks_stats_zi <- ks_stats[191:380]
ks_stats_risk <- ks_stats[1:190]

summary(ks_stats_risk)
summary(ks_stats_zi)

histbreaks <- seq(0,.24,by=.02)

DEFAULTTHEME <- theme_classic() +
  theme(text = element_text(size = 28))

zi_ks_plot <- ggplot(data.frame(x=ks_stats_zi),aes(x=x)) +
  DEFAULTTHEME +
  geom_histogram(aes(y=..density..),breaks = histbreaks,fill='transparent',colour='black') +
  labs(x='',y='') +
  coord_cartesian(xlim=c(0,.2))
ggsave(plot=zi_ks_plot,file = file.path(plotpath,paste0("loaloa-ks-zi.pdf")),width=7,height=7)

risk_ks_plot <- ggplot(data.frame(x=ks_stats_risk),aes(x=x)) +
  DEFAULTTHEME +
  geom_histogram(aes(y=..density..),breaks = histbreaks,fill='transparent',colour='black') +
  labs(x='',y='') +
  coord_cartesian(xlim=c(0,.2))
ggsave(plot=risk_ks_plot,file = file.path(plotpath,paste0("loaloa-ks-risk.pdf")),width=7,height=7)

cat(paste0("All done! Go to: ",plotpath," to see the results.\n"))