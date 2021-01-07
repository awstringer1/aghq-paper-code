### SIR model using AGHQ ###

### Epidemic Modelling with epiLMCT ###
# install.packages("EpiILMCT")
library(tidyverse)
library(TMB)
library(EpiILMCT)
library(plot3D)
library(trustOptim)
library(Matrix)
library(parallel)
options(mc.cores = parallel::detectCores())
library(aghq)

globalpath <- "~/phd/projects/best-friends-gang/normalizing-constant/"
plotpath <- paste0(globalpath,"figures/disease/")

# Get the Tomato disease data
data("tswv", package = "EpiILMCT")

## TMB ##

precompile()
setwd("~/phd/projects/best-friends-gang/normalizing-constant/code/")

# Compile-- only need to do if source has changed
# compile("02_disease.cpp") # Optimized
# compile("02_disease.cpp","-O0 -g") # Unoptimized, accurate line numbers for debugging

dyn.load(dynlib("02_disease"))

# Create the functions
dat <- tswv$tswvsir
dat$epidat <- dat$epidat[order(dat$epidat[ ,4]), ]

I <- dat$epidat[ ,4]
R <- dat$epidat[ ,2]
infected <- !is.infinite(I)

datlist <- list(
  D = as.matrix(dist(dat$location[dat$epidat[ ,1], ])), 
  I = I,
  R = R,
  infected = as.numeric(infected[infected])
)

ff <- MakeADFun(data = datlist,
                parameters = list(theta1 = 0,theta2 = 0),
                DLL = "02_disease",
                ADreport = FALSE,
                silent = TRUE)


# AGHQ ----

resultslist <- list(
  aghq(ff,3,c(0,0)),
  aghq(ff,5,c(0,0)),
  aghq(ff,7,c(0,0)),
  aghq(ff,9,c(0,0)),
  aghq(ff,11,c(0,0)),
  aghq(ff,13,c(0,0))
)
names(resultslist) <- seq(3,13,by=2)

# Timing

microbenchmark::microbenchmark(
  aghq(ff,3,c(0,0)),
  aghq(ff,5,c(0,0)),
  aghq(ff,7,c(0,0)),
  aghq(ff,9,c(0,0)),
  aghq(ff,11,c(0,0)),
  aghq(ff,13,c(0,0))
)

# Unit: milliseconds
# expr       min        lq      mean    median        uq       max neval     cld
# aghq(ff, 3, c(0, 0))  212.8668  213.8543  215.2693  214.5139  215.2885  224.8532   100 a      
# aghq(ff, 5, c(0, 0))  272.1243  273.3979  275.2554  274.2651  275.2142  285.1065   100  b     
# aghq(ff, 7, c(0, 0))  360.7452  362.9744  365.2792  364.2050  365.9398  374.8743   100   c    
# aghq(ff, 9, c(0, 0))  479.5925  482.7099  486.3982  484.4351  485.8624  662.2081   100    d   
# aghq(ff, 11, c(0, 0))  627.7888  632.7916  635.9133  634.4418  637.3114  650.6599   100     e  
# aghq(ff, 13, c(0, 0))  808.3682  813.1594  816.0081  815.2659  818.2919  830.6486   100      f 

.215 * 110000/3408
.275 * 110000/3408
.365 * 110000/3408
.486 * 110000/3408
.636 * 110000/3408
.815 * 110000/3408

# Get the MCMC results
# nsim <- 110000
# nchains <- 4
# 
# covsus <- list(NULL)
# covsus[[1]] <- list(0.02, c("gamma", 1, 0.01, 0.01))
# covsus[[2]] <- rep(1, length(tswv$tswvsir$epidat[,1]))
# kernel1 <- list(2, c("gamma", 1, 0.01, 0.1))
# set.seed(524837)
# system.time({
#   tswv.full.observed <- epictmcmc(object = tswv$tswvsir, 
#                                   distancekernel = "powerlaw",
#                                   datatype = "known epidemic",
#                                   nsim = nsim, 
#                                   nchains = nchains,
#                                   parallel = TRUE,
#                                   control.sus = covsus, 
#                                   kernel.par = kernel1)
# })
# plot(tswv.full.observed,plottype = "parameter", start = 10001, thin = 10,density = FALSE)
# # 110000 iter, 4 chains in parallel: 3408 secs (about an hour)
# save(tswv.full.observed,file = paste0(globalpath,"results/sir-mcmc-20200717.RData"))
load(paste0(globalpath,"results/sir-mcmc-20200717.RData"))

postsamps <- tswv.full.observed$parameter.samples %>%
  map(as.data.frame) %>%
  map(as_tibble) %>%
  map(~slice(.x,10001:110000)) %>%
  reduce(bind_rows) %>%
  rename(alpha = `Alpha_s[1]`,beta = `Spatial parameter`)

postmeans <- apply(postsamps,2,mean)
postsd <- apply(postsamps,2,sd)
postquantiles2.5 <- apply(postsamps,2,function(x) quantile(x,.025))
postquantiles97.5 <- apply(postsamps,2,function(x) quantile(x,.975))

mcmcplotalpha <- postsamps %>%
  ggplot(aes(x = alpha)) +
  theme_classic() +
  geom_histogram(aes(y = ..density..),bins = 100,colour = "transparent",fill = "grey",alpha = .5)

mcmcplotbeta <- postsamps %>%
  ggplot(aes(x = beta)) +
  theme_classic() +
  geom_histogram(aes(y = ..density..),bins = 100,colour = "transparent",fill = "grey",alpha = .5)

get_plots <- function(K) {
  # Get:
  # plots for alpha, beta
  # moments and quantiles
  # KS
  results <- resultslist[[as.character(K)]]
  
  plotpointsalpha <- results$marginals[[1]]
  plotpointsbeta <- results$marginals[[2]]
  
  plotdatalpha <- compute_pdf_and_cdf(plotpointsalpha,list(totheta = log,fromtheta = exp))
  plotdatbeta <- compute_pdf_and_cdf(plotpointsbeta,list(totheta = log,fromtheta = exp))
  
  alphaplot <- mcmcplotalpha +
    geom_line(
      data = plotdatalpha,
      mapping = aes(x = transparam,y = pdf_transparam)
    ) +
    geom_point(
      data = plotpointsalpha,
      mapping = aes(x = exp(theta1),y = exp(logmargpost - theta1)), # Jacobian
      size = 5
    ) +
    labs(title = "",x = "",y = "Density") +
    theme(text = element_text(size = 28)) +
    coord_cartesian(xlim = c(0.005,0.0255),ylim = c(0,200))
  
  betaplot <- mcmcplotbeta +
    geom_line(
      data = plotdatbeta,
      mapping = aes(x = transparam,y = pdf_transparam)
    ) +
    geom_point(
      data = plotpointsbeta,
      mapping = aes(x = exp(theta2),y = exp(logmargpost - theta2)), # Jacobian
      size = 5
    ) +
    labs(title = "",x = "",y = "Density") +
    theme(text = element_text(size = 28)) +
    coord_cartesian(xlim = c(0.8,2.05),ylim = c(0,3))
  
  
  list(alpha = alphaplot,beta = betaplot)
}

theplots <- list(
  get_plots(3),
  get_plots(5),
  get_plots(7)
)
names(theplots) <- seq(3,7,by=2)

# Save plots to disk

ggsave(paste0(plotpath,"disease-alpha-K3-20201012.pdf"),
       theplots[["3"]]$alpha,
       width = 7,height = 7)
ggsave(paste0(plotpath,"disease-alpha-K5-20201012.pdf"),
       theplots[["5"]]$alpha,
       width = 7,height = 7)
ggsave(paste0(plotpath,"disease-alpha-K7-20201012.pdf"),
       theplots[["7"]]$alpha,
       width = 7,height = 7)

ggsave(paste0(plotpath,"disease-beta-K3-20201012.pdf"),
       summaries[["3"]]$beta,
       width = 7,height = 7)
ggsave(paste0(plotpath,"disease-beta-K5-20201012.pdf"),
       summaries[["5"]]$beta,
       width = 7,height = 7)
ggsave(paste0(plotpath,"disease-beta-K7-20201012.pdf"),
       summaries[["7"]]$beta,
       width = 7,height = 7)

# Summaries
get_summaries <- function(K) {
  results <- resultslist[[as.character(K)]]
  
  themeans <- compute_moment(results$normalized_posterior,exp)
  thesds <- compute_moment(results$normalized_posterior,function(x) (exp(x) - themeans)^2) %>% sqrt()
  alphaquants <- exp(compute_quantiles(results$marginals[[1]]))
  betaquants <- exp(compute_quantiles(results$marginals[[2]]))
  
  tibble(
    alphamean = themeans[1],
    betamean = themeans[2],
    alphasd = thesds[1],
    betasd = thesds[2],
    alpha2.5 = alphaquants[1],
    beta2.5 = betaquants[1],
    alpha97.5 = alphaquants[2],
    beta97.5 = betaquants[2]
  )
}

map(seq(3,13,by=2),get_summaries) %>% reduce(bind_rows)


