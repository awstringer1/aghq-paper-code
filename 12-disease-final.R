### SIR model using AGHQ ###

### Epidemic Modelling with epiLMCT ###
# install.packages("EpiILMCT")
library(tidyverse)
library(TMB)
precompile()
library(tmbstan)
library(Matrix)
library(parallel)
options(mc.cores = parallel::detectCores())
library(aghq)

set.seed(573489)

globalpath <- normalizePath(tempdir(),winslash='/')
plotpath <- normalizePath(file.path(globalpath,"disease"),winslash='/',mustWork = FALSE)
if (!dir.exists(plotpath)) dir.create(plotpath)
plotstamp <- '-2021-05-12'

# Get the Tomato disease data
data("tswv", package = "EpiILMCT")

## TMB ##
# get the template from the aghq package
file.copy(normalizePath(system.file('extsrc/02_disease.cpp',package='aghq'),winslash='/'),globalpath)
compile(normalizePath(file.path(globalpath,'02_disease.cpp'),winslash='/'))
dyn.load(normalizePath(dynlib(file.path(globalpath,"02_disease")),winslash='/'))


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

cntrl <- default_control(negate = TRUE,method_summaries = 'correct')
trans <- make_transformation(log,exp)
resultslist <- list(
  aghq(ff,3,c(0,0),control = cntrl,transformation = trans),
  aghq(ff,5,c(0,0),control = cntrl,transformation = trans),
  aghq(ff,7,c(0,0),control = cntrl,transformation = trans),
  aghq(ff,9,c(0,0),control = cntrl,transformation = trans),
  aghq(ff,11,c(0,0),control = cntrl,transformation = trans),
  aghq(ff,13,c(0,0),control = cntrl,transformation = trans)
)
names(resultslist) <- seq(3,13,by=2)

# Timing

microbenchmark::microbenchmark(
  aghq(ff,3,c(0,0),control = cntrl),
  aghq(ff,5,c(0,0),control = cntrl),
  aghq(ff,7,c(0,0),control = cntrl),
  aghq(ff,9,c(0,0),control = cntrl),
  aghq(ff,11,c(0,0),control = cntrl),
  aghq(ff,13,c(0,0),control = cntrl)
)

# Unit: milliseconds
# expr                                      min       lq     mean   median       uq      max neval
# aghq(ff, 3, c(0, 0), control = cntrl) 100.6553 100.8632 101.1796 100.9300 101.0215 105.3452   100
# aghq(ff, 5, c(0, 0), control = cntrl) 159.8424 160.0872 160.4006 160.2457 160.4171 164.9010   100
# aghq(ff, 7, c(0, 0), control = cntrl) 223.3692 223.8359 224.6426 223.9588 224.2837 230.4604   100
# aghq(ff, 9, c(0, 0), control = cntrl) 284.8110 285.2891 286.3028 285.4271 285.8088 297.2387   100
# aghq(ff, 11, c(0, 0), control = cntrl) 356.8188 357.1761 358.2884 357.3310 358.6568 364.1082   100
# aghq(ff, 13, c(0, 0), control = cntrl) 441.1805 441.6348 442.7299 441.7533 442.1720 463.5870   100

# Get the MCMC results
numiter <- 1e04
numwarmup <- 1e03
stanmod <- tmbstan(
  ff,
  chains = 4,
  cores = 4,
  iter = numiter,
  warmup = numwarmup,
  init = c(0,0),
  seed = 124698,
  algorithm = "NUTS"
)

# Timings
mcmctime <- max(apply(rstan::get_elapsed_time(stanmod),1,sum))
round(c(
  .101 * numiter/mcmctime,
  .160 * numiter/mcmctime,
  .224 * numiter/mcmctime,
  .285 * numiter/mcmctime,
  .357 * numiter/mcmctime,
  .441 * numiter/mcmctime
),0)


postsamps <- as.data.frame(stanmod)
postsamps$alpha <- exp(postsamps[ ,1])
postsamps$beta <- exp(postsamps[ ,2])
postsamps <- postsamps[ ,c('alpha','beta')]

postmeans <- apply(postsamps,2,mean)
postsd <- apply(postsamps,2,sd)
postquantiles2.5 <- apply(postsamps,2,function(x) quantile(x,.025))
postquantiles97.5 <- apply(postsamps,2,function(x) quantile(x,.975))

mcmcplotalpha <- postsamps %>%
  ggplot(aes(x = alpha)) +
  theme_classic() +
  geom_histogram(aes(y = ..density..),bins = 50,colour = "transparent",fill = "grey",alpha = .5)

mcmcplotbeta <- postsamps %>%
  ggplot(aes(x = beta)) +
  theme_classic() +
  geom_histogram(aes(y = ..density..),bins = 50,colour = "transparent",fill = "grey",alpha = .5)

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


ggsave(file.path(plotpath,paste0("disease-alpha-K3-",plotstamp,".pdf")),
       theplots[["3"]]$alpha,
       width = 7,height = 7)
ggsave(file.path(plotpath,paste0("disease-alpha-K5-",plotstamp,".pdf")),
       theplots[["5"]]$alpha,
       width = 7,height = 7)
ggsave(file.path(plotpath,paste0("disease-alpha-K7-",plotstamp,".pdf")),
       theplots[["7"]]$alpha,
       width = 7,height = 7)

ggsave(file.path(plotpath,paste0("disease-beta-K3-",plotstamp,".pdf")),
       theplots[["3"]]$beta,
       width = 7,height = 7)
ggsave(file.path(plotpath,paste0("disease-beta-K5-",plotstamp,".pdf")),
       theplots[["5"]]$beta,
       width = 7,height = 7)
ggsave(file.path(plotpath,paste0("disease-beta-K7-",plotstamp,".pdf")),
       theplots[["7"]]$beta,
       width = 7,height = 7)

# Summaries
get_summaries <- function(K) {
  results <- resultslist[[as.character(K)]]

  themeans <- compute_moment(results,nn=1,method='correct')
  thesds <- sqrt(compute_moment(results,nn=2,method='correct',type='central'))

  quants <- compute_quantiles(results)
  alphaquants <- quants[[1]]
  betaquants <- quants[[2]]

  aghqpostsamps <- sample_marginal(results,nrow(postsamps),transformation = default_transformation())

  suppressWarnings({
    ks_alpha <- ks.test(
      postsamps$alpha,
      exp(aghqpostsamps[[1]])
    )$statistic

    ks_beta <- ks.test(
      postsamps$beta,
      exp(aghqpostsamps[[2]])
    )$statistic
  })

  tibble(
    alphamean = themeans[1],
    betamean = themeans[2],
    alphasd = thesds[1],
    betasd = thesds[2],
    alpha2.5 = alphaquants[1],
    beta2.5 = betaquants[1],
    alpha97.5 = alphaquants[2],
    beta97.5 = betaquants[2],
    ks_alpha = ks_alpha,
    ks_beta = ks_beta
  )
}

postsummaries <- map(seq(3,13,by=2),get_summaries) %>% reduce(bind_rows)
readr::write_csv(postsummaries,file.path(plotpath,paste0("posteriorsummarytable-",plotstamp,".csv")))
knitr::kable(
  postsummaries,
  digits = 5,
  format = 'markdown'
)




