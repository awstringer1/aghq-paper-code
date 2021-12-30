### Loa Loa Zip Model ###

## Load the data ##
library(TMB)
precompile()
library(geostatsp)
data(loaloa,package = "geostatsp")
library(aghq)
library(tmbstan)

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

savestamp <- "20211230-v1"
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

if (doaghq) {
  tm <- Sys.time()
  cat("Doing AGHQ, time = ",format(tm),"\n")
  loaloazipquad <- aghq::marginal_laplace_tmb(
    ff,
    3,
    startingvalue = c(paraminit$logkappa,paraminit$logtau)
  )
  aghqtime <- difftime(Sys.time(),tm,units = 'secs')
  save(loaloazipquad,file = file.path(savepath,paste0("loaloazipquad",savestamp,".RData")))
  cat("AGHQ took: ",format(aghqtime),"\n")
}

## MCMC ----
if (domcmc){
  tm <- Sys.time()
  cat("Doing MCMC, time = ",format(tm),"\n")
  loaloazipmcmc <- tmbstan(
    ff,
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
  
  # NOTE: for run the week of 15/04/2021, got the following warnings with default
  # settings and 8 chains with 10,000 iterations (incl 1,000 warmups):
  # There were 238 divergent transitions after warmup.
  # There were 673 transitions after warmup that exceeded the maximum treedepth.
  pdf(file = paste0(plotpath,"mcmc-pairs-plot",savestamp,".pdf"))
  pairs(loaloazipmcmc,pars = c("W[191]","W[381]","logkappa","logtau"))
  dev.off()
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
  postU <- loazippostsamples$samps[c(1:190), ]
  postV <- loazippostsamples$samps[c(192:381), ]
  
  postBeta <- loazippostsamples$samps[c(191,382), ]
  
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
  
  fieldbrickzip_withint <- fieldbrickzip + postBeta[1, ]
  fieldbrickrisk_withint <- fieldbrickrisk + postBeta[2, ]
  
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
  # Only do 100
  idx <- sample.int(nrow(samps),100)
  
  postU <- t(samps[idx,c(1:190)])
  postV <- t(samps[idx,c(192:381)])
  
  postBeta <- samps[idx,c(191,382)]
  
  theta <- samps[idx,c(383,384)]
  colnames(theta) <- c('logkappa','logtau')
  
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
  
  fieldbrickzip_withint <- fieldbrickzip + postBeta[ ,1]
  fieldbrickrisk_withint <- fieldbrickrisk + postBeta[ ,2]
  
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

# Do the pairs plots
if (dopostsamplingaghq & dopostsamplingmcmc) {
  set.seed(5678798)
  # MCMC
  mcmcbetasamps <- as.data.frame(loaloazipmcmc)[ ,c(191,382)]
  pdf(file = file.path(plotpath,paste0("mcmcpairsplot-",savestamp,".pdf")),width = 7,height = 7)
  pairs(loaloazipmcmc,pars = c("W[191]","W[382]"),text.panel = function(x,y,labels,cex,font,...) NULL)
  dev.off()
  
  
  
  # AGHQ
  aghqbetasamps <- t(sample_marginal(loaloazipquad,nrow(mcmcbetasamps))$samps[c(191,382), ])
  pdf(file = file.path(plotpath,paste0("aghqpairsplot-",savestamp,".pdf")),width = 7,height = 7)
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
  }
  smoothscatterforhere <- function(x,y,...) {
    smoothScatter(x,y,add = TRUE,nrpoints = 0)
  }
  pairs(aghqbetasamps,
        panel = smoothscatterforhere,
        diag.panel = panel.hist,
        text.panel = function(x,y,labels,cex,font,...) NULL)
  dev.off()
  
  
  
}

# Write the timing table
timingtable <- data.frame(
  task = c("AGHQ","MCMC","AGHQSim","MCMCSim"),
  time = c(aghqtime,mcmctime,aghqsimtime,mcmcsimtime)
)
readr::write_csv(timingtable,file.path(savepath,paste0("timing-table",savestamp,".csv")))
cat("All done!\n")