### Loa Loa zero-inflated example for AGHQ paper ###

## Load the data ##
library(tidyverse)
library(TMB)
library(geostatsp)
library(aghq)
data(loaloa,package = "geostatsp")

# Set working dir for TMB executation
setwd("/storage/phd/projects/best-friends-gang/normalizing-constant/code/")

# Create covariate rasters
elevationLoa <- elevationLoa - 750
elevLow <- reclassify(elevationLoa, c(0, Inf, 0))
elevHigh <- reclassify(elevationLoa, c(-Inf, 0, 0))

rcl <- rbind(c(9, 8), c(5, 2), c(11, 2), c(12, 14), c(13, 14))
ltLoaRe <- reclassify(ltLoa, rcl)
levels(ltLoaRe) = levels(ltLoa)

## Prepare the "inner" model ##

model_data <- new.env()

# Design matrices
model_data$Amat <- Diagonal(nrow(loaloa))

model_data$Xmat <- NULL
# Design matrix: zip model and risk model are the same
model_data$design <- bdiag(
  # ZIP
  cbind(
    model_data$Amat
  ),
  # Risk
  cbind(
    model_data$Amat  
  )
)

# Response
model_data$y <- loaloa@data$y
model_data$N <- loaloa@data$N

## Dimensions
model_data$n <- nrow(model_data$Xmat) # Number of obs
model_data$p <- ncol(model_data$Xmat) # Number of betas
model_data$m <- ncol(model_data$Amat) * 2 # Number of spatial points
model_data$Wd <- ncol(model_data$design) # Number of total params
# Check
stopifnot(model_data$Wd == model_data$m + model_data$p)
# FOUR matern parameters, two spatial processes
model_data$S <- 2

## Prior distributions ##
# Use the same prior for both sets of Matern params
sigma_u <- 4
# sigma_u <- 1
sigma_alpha <- .025
rho_u <- 5e05
# rho_u <- 5e04
rho_alpha <- .975

# PC Prior for kappa,tau
model_data$maternconstants$d <- 2 # Dimension of spatial field, fixed
model_data$maternconstants$nu <- 2 # Shape param, fixed
get_kappa <- function(sigma,rho) sqrt(8*model_data$maternconstants$nu)/rho
get_tau <- function(sigma,rho) sigma * get_kappa(sigma,rho)^(model_data$maternconstants$nu) * sqrt(gamma(model_data$maternconstants$nu + model_data$maternconstants$d/2) * (4*pi)^(model_data$maternconstants$d/2) / gamma(model_data$maternconstants$nu))
get_sigma <- function(kappa,tau) tau / (kappa^(model_data$maternconstants$nu) * sqrt(gamma(model_data$maternconstants$nu + model_data$maternconstants$d/2) * (4*pi)^(model_data$maternconstants$d/2) / gamma(model_data$maternconstants$nu)))
get_rho <- function(kappa,tau) sqrt(8*model_data$maternconstants$nu) / kappa

# Log prior for ONE set of theta = (log(kappa),log(tau))
log_prior_theta <- function(theta) {
  # theta = (log(kappa),log(tau))
  kappa <- exp(theta[1])
  tau <- exp(theta[2])
  
  lambda1 <- -(rho_u / sqrt(8*model_data$maternconstants$nu)) ^ (model_data$maternconstants$d/2) * log(rho_alpha)
  lambda2 <- -kappa^(-model_data$maternconstants$nu) * sqrt( gamma(model_data$maternconstants$nu) / ( gamma(model_data$maternconstants$nu + model_data$maternconstants$d/2) * (4*pi)^(model_data$maternconstants$d/2) ) ) * log(sigma_alpha) / sigma_u
  
  log(model_data$maternconstants$d) - log(2) + log(lambda1) + log(lambda2) + (model_data$maternconstants$d/2 - 1) * log(kappa) - lambda1 * kappa^(model_data$maternconstants$d/2) - lambda2 * tau + sum(theta)
}

model_data$beta_prec <- .001

Q_matrix <- function(theta) {
  # theta = log(kappa), log(tau)
  
  kappa <- as.numeric(unname(exp(theta[1])))
  tau <- as.numeric(unname(exp(theta[2])))
  
  sig <- get_sigma(kappa,tau)
  rho <- get_rho(kappa,tau)
  
  # Matern
  mm <- geostatsp::matern(
    loaloa,
    param = c("variance" = sig^2,"range" = rho,"shape" = 2),
    type = "precision"
  )
  
  
  bdiag(mm,mm)
}

log_prior_W <- function(W,theta,Q = NULL) {
  if (is.null(Q)) Q <- Q_matrix(theta)
  -(1/2) * as.numeric(crossprod(W,crossprod(Q,W))) + .5 * determinant(Q,logarithm=TRUE)$modulus
}


## Likelihood ----

# compile("10_loaloazip.cpp")
dyn.load(dynlib("10_loaloazip"))

datlist <- list(
  y = model_data$y,
  N = model_data$N
)
paraminit <- list(
  eta1 = rep(0,nrow(loaloa)),
  eta2 = rep(0,nrow(loaloa))
)

ll <- MakeADFun(data = datlist,
                parameters = paraminit,
                DLL = "10_loaloazip",
                ADreport = FALSE,
                silent = TRUE)

# Verified likelihood and derivatives are correct...




## Posterior

make_eta <- function(W) as.numeric(model_data$design %*% W)

log_posterior_W <- function(W,theta,Q = NULL) {
  eta <- make_eta(W)
  if (is.null(Q)) Q <- Q_matrix(theta)
  log_prior_W(W,theta,Q) + ll$fn(eta)
}

grad_log_posterior_W <- function(W,theta,Q = NULL) {
  eta <- make_eta(W)
  if (is.null(Q)) Q <- Q_matrix(theta)
  as.numeric(-Q %*% cbind(W) + t(model_data$design) %*% t(ll$gr(eta)))
}

H_matrix <- function(W,theta,Q = NULL) {
  # minus the hessian of the log posterior
  eta <- make_eta(W)
  CE <- -1 * ll$he(eta)
  if (is.null(Q)) Q <- Q_matrix(theta)
  
  CW <- crossprod(model_data$design,crossprod(CE,model_data$design))
  
  as(Q + CW,"dgCMatrix")
}

## Quadrature ##

log_posterior_joint <- function(W,theta) {
  # theta1 <- theta[1:2]
  # theta2 <- theta[3:4]
  theta1 <- theta
  log_prior_theta(theta1) +
    # log_prior_theta(theta2) +
    log_posterior_W(W,theta)
}

ff <- list(
  fn = log_posterior_joint,
  gr = grad_log_posterior_W,
  he = function(W,theta) -1 * H_matrix(W,theta)
)

startingsig <- .988
startingrho <- 4.22*1e04
startingtheta <- rep(c(
  log(get_kappa(startingsig,startingrho)),
  log(get_tau(startingsig,startingrho))
),2)
startingtheta <- startingtheta[1:2]

# Get starting values for W
ff2 <- list(
  fn = function(W) ff$fn(W,startingtheta),
  gr = function(W) ff$gr(W,startingtheta),
  he = function(W) ff$he(W,startingtheta)
)

opt <- aghq::optimize_theta(ff2,rep(0,model_data$Wd),control = list(method = 'sparse_trust'))

tm <- Sys.time()
loaloazipquad <- aghq::marginal_laplace(ff,3,list(W = opt$mode,theta = startingtheta),control = list(method = 'trust',inner_method = 'sparse_trust'))
runtime <- difftime(Sys.time(),tm,units='secs')
cat("Run time for loaloa AGHQ:",runtime,"seconds.\n")
save(loaloazipquad,file = "loaloazipquad-20210208-constrequal-prevpriors.RData")
# load(file = "loaloazipquad-20210208-constrequal-prevpriors.RData")


## Sampling ##
simulate_spatial_fields <- function(U,theta,pointsdata,resolution = list(nrow = 100,ncol = 100)) {
  # U: matrix of samples, each column is a sample
  # theta: tibble of theta values
  # Draw from U*|U
  fieldlist <- vector(mode = 'list',length = nrow(theta))
  for (i in 1:length(fieldlist)) {
    fielddat <- pointsdata
    fielddat@data <- data.frame(w = as.numeric(U[ ,i]))
    
    # Back-transform the Matern params
    kappa <- exp(theta$theta1[i])
    tau <- exp(theta$theta2[i])
    sig <- get_sigma(kappa,tau)
    rho <- get_rho(kappa,tau)
    # Simulate from the two fields
    capture.output({
      fieldlist[[i]] <- geostatsp::RFsimulate(
        model = c("variance" = sig^2,"range" = rho,"shape" = 2),
        data = fielddat,
        x = raster(fielddat,nrows = resolution$nrow,ncols = resolution$ncol),
        n = 1
      )
    })
  }
  brick(fieldlist)
}


# Post samples
set.seed(80979)
loazippostsamples <- aghq::sample_marginal(loaloazipquad,100)
nn <- nrow(loazippostsamples$samps)
UUzip <- loazippostsamples$samps[1:(nn/2), ]
UUrisk <- loazippostsamples$samps[(1+nn/2):nn, ]

tm <- Sys.time()
fieldbrickzip <- simulate_spatial_fields(
  UUzip,
  loazippostsamples$theta,
  loaloa,
  resolution = list(nrow = 200,ncol = 400)
)
runtime <- difftime(Sys.time(),tm,units='secs')
cat("Run time for fieldbrickzip:",runtime,"seconds.\n")

tm <- Sys.time()
fieldbrickrisk <- simulate_spatial_fields(
  UUrisk,
  loazippostsamples$theta,
  loaloa,
  resolution = list(nrow = 200,ncol = 400)
)
runtime <- difftime(Sys.time(),tm,units='secs')
cat("Run time for fieldbrickzip:",runtime,"seconds.\n")

save(loazippostsamples,file = "loaloazipquad-20210208-postsamples.RData")
raster::writeRaster(fieldbrickzip,file = "loaloazipquad-20210208-fieldbrickzip.grd",overwrite = FALSE)
raster::writeRaster(fieldbrickrisk,file = "loaloazipquad-20210208-fieldbrickrisk.grd",overwrite = FALSE)

# load(file = "loaloazipquad-20210208-postsamples.RData")
# fieldbrickzip <- raster::brick("loaloazipquad-20210208-fieldbrickzip.grd")
# fieldbrickrisk <- raster::brick("loaloazipquad-20210208-fieldbrickrisk.grd")


# Posterior means
simfieldsmeanzip <- calc(fieldbrickzip,function(x) mean(1 / (1 + exp(-x))))
simfieldsmeanrisk <- calc( (1/(1+exp(-fieldbrickzip))) * (1 /(1+exp(-fieldbrickrisk))),mean)

# Exceedence probabilities


# Border for Cameroon and Nigeria
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

# Mean of U
plotraster <- simfieldsmeanzip

predcols <- mapmisc::colourScale(
  plotraster,
  breaks = c(0,.05,.1,.5,.6,.7,.8,.9,.95,1),
  style = "fixed",
  col = "Spectral",
  rev = TRUE
  # dec = -log10(0.05)
)

plotraster <- mask(plotraster,fullborderouter)

# par(mfrow = c(2,1))

pdf(file = 'loaloazip-suitprob.pdf')
mapmisc::map.new(loaloa)
plot(plotraster,
     col = predcols$col,
     breaks = predcols$breaks,
     legend=FALSE, add=TRUE)
points(loaloa,pch = '.')
plot(fullborder,add = TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
plot(fullborderouter,add = TRUE)
mapmisc::legendBreaks('right', predcols, cex=1, bty='o',inset = 0.01)
dev.off()

# Mean of lam
plotraster <- simfieldsmeanrisk

predcols <- mapmisc::colourScale(
  plotraster,
  breaks = c(0,.1,.15,.2,.25,.3,.4,.5,.6),
  style = "fixed",
  col = "Spectral",
  rev = TRUE
  # dec = -log10(0.05)
)

plotraster <- mask(plotraster,fullborderouter)

pdf(file = 'loaloazip-infectprob.pdf')
mapmisc::map.new(plotraster)
plot(plotraster,
     col = predcols$col,
     breaks = predcols$breaks,
     legend=FALSE, add=TRUE)
points(loaloa,pch = '.')
plot(fullborder,add = TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
plot(fullborderouter,add = TRUE)
mapmisc::legendBreaks('right', predcols, cex=1, bty='o',inset = .01)
dev.off()
