### Simulations for AGHQ paper ###
# Note: these simulations predate the aghq package, so this package is not used here.

## Setup ----
library(tidyverse)
library(trustOptim)
library(Matrix)
library(parallel)
options(mc.cores = parallel::detectCores())

set.seed(8907342)

# Global constants
# globalpath <- "~/phd/projects/best-friends-gang/normalizing-constant/"
globalpath <- normalizePath(tempdir(),winslash='/')
plotpath <- normalizePath(file.path(globalpath,"figures"),winslash='/')
if (!dir.exists(plotpath)) dir.create(plotpath)

# Function to approximate

logft <- function(lambda,y) sum(y) * log(lambda) - (length(y) + 1) * lambda - sum(lgamma(y+1))

etafromlambda <- function(lambda) log(lambda)
lambdafrometa <- function(eta) exp(eta)
detadlambda <- function(lambda) 1/lambda
dlambdadeta <- function(eta) exp(eta)

logfteta <- function(eta,y) {
  lambda <- lambdafrometa(eta)
  logft(lambda,y) + log(abs(dlambdadeta(eta)))
}

truelogint <- function(y) lgamma(1 + sum(y)) - (1 + sum(y)) * log(length(y) + 1) - sum(lgamma(y+1))

aghq <- function(K,y) {
  n <- length(y)
  # This one has a closed form answer
  eta_hat <- etafromlambda((sum(y) + 1 )/ (length(y)+1))
  
  thehess <- -1 * numDeriv::hessian(function(eta) logfteta(eta,y),eta_hat)
  
  intgrid <- mvQuad::createNIGrid(1,"GHe",K)
  mvQuad::rescale(intgrid,m = eta_hat,C = solve(thehess))
  
  # Do logsumexp and return log of int const
  # Relative error will then be on the log scale
  nn <- as.numeric(mvQuad::getNodes(intgrid))
  ll <- logfteta(nn,y)
  ww <- as.numeric(mvQuad::getWeights(intgrid))
  logintconst <- matrixStats::logSumExp(ll + log(ww))
  
  list(
    logintconst = logintconst
  )
}

## Simulations ----

ntodo <- seq(1,100,by=1)

ktodo <- c(3,5,7,11)

sims <- expand.grid(
  ntodo,
  ktodo
)
sims <- as_tibble(sims)
colnames(sims) <- c("n","k")
sims$logrelerror <- 0
sims$rate <- 0
sims$diffoflogs <- 0
simlist <- split(sims,1:nrow(sims))

do_sim <- function(lst) {
  n <- lst$n
  k <- lst$k
  lambda <- 5 # Set these here, for running in parallel
  # M <- 100
  M <- 100
  out <- slice(lst,rep(1,M))
  for (i in 1:M) {
    y <- rpois(n,lambda)
    logintconst <- aghq(k,y)$logintconst
    truelogintconst <- truelogint(y)
    out[i,"logrelerror"] <- log(abs(exp(truelogintconst - logintconst) - 1))
    rk <- floor((k+2)/3)
    out[i,"rate"] <- rk
    out[i,"diffoflogs"] <- out[i,"logrelerror"] + rk * log(n)
  }
  out
}

thesims <- mclapply(simlist,do_sim)
thesimframe <- reduce(thesims,bind_rows)


simplotk3 <- thesimframe %>%
  filter(k == 3) %>%
  ggplot(aes(x = n,y = diffoflogs)) +
  theme_classic() +
  geom_point(pch = 21) +
  scale_x_continuous(breaks = seq(0,100,by=10)) +
  theme(text = element_text(size = 28)) +
  labs(y = "")

simplotk5 <- thesimframe %>%
  filter(k == 5) %>%
  ggplot(aes(x = n,y = diffoflogs)) +
  theme_classic() +
  geom_point(pch = 21) +
  scale_x_continuous(breaks = seq(0,100,by=10)) +
  theme(text = element_text(size = 28)) +
  labs(y = "")

simplotk7 <- thesimframe %>%
  filter(k == 7) %>%
  ggplot(aes(x = n,y = diffoflogs)) +
  theme_classic() +
  geom_point(pch = 21) +
  scale_x_continuous(breaks = seq(0,100,by=10)) +
  theme(text = element_text(size = 28)) +
  labs(y = "")

simplotk11 <- thesimframe %>%
  filter(k == 11) %>%
  ggplot(aes(x = n,y = diffoflogs)) +
  theme_classic() +
  geom_point(pch = 21) +
  scale_x_continuous(breaks = seq(0,100,by=10)) +
  theme(text = element_text(size = 28)) +
  labs(y = "")

ggsave(
  filename = file.path(plotpath,paste0("relratescatterplot3.pdf")),
  plot = simplotk3,
  width = 7,
  height = 7)

ggsave(
  filename = file.path(plotpath,paste0("relratescatterplot5.pdf")),
  plot = simplotk5,
  width = 7,
  height = 7)

ggsave(
  filename = file.path(plotpath,paste0("relratescatterplot7.pdf")),
  plot = simplotk7,
  width = 7,
  height = 7)

ggsave(
  filename = file.path(plotpath,paste0("relratescatterplot11.pdf")),
  plot = simplotk11,
  width = 7,
  height = 7)