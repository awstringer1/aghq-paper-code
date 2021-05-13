### Milky Way mass estimation using AGHQ ###

library(tidyverse)
library(patchwork)
library(TMB)
library(tmbstan)
library(ipoptr)
library(aghq)
library(Matrix)

data(gcdatalist,package = 'aghq')
precompile()

savestamp <- "20210417-v1"

set.seed(4365789)

globalpath <- tempdir()
plotpath <- file.path(globalpath,"astro")
if (!dir.exists(plotpath)) dir.create(plotpath)
plotstamp <- '-2021-05-12'

# Get the Tomato disease data
data("tswv", package = "EpiILMCT")

## TMB ##
# get the template from the aghq package
file.copy(system.file('extsrc/01_astro.cpp',package='aghq'),globalpath)
compile(file.path(globalpath,'01_astro.cpp'))
dyn.load(dynlib(file.path(globalpath,"01_astro")))

## Parameter transformations ##
parambounds <- list(
  Psi0 = c(1,200),
  gamma = c(.3,.7),
  alpha = c(3.0,3.7),
  beta = c(-.5,1)
)

get_psi0 <- function(theta) {
  # theta = -log( (Psi0 - 1) / (200 - 1) )
  (parambounds$Psi0[2] - parambounds$Psi0[1]) * 
    exp(-exp(theta)) + parambounds$Psi0[1]
}
get_theta1 <- function(Psi0) log(
  -log( 
    (Psi0 - parambounds$Psi0[1]) / (parambounds$Psi0[2] - parambounds$Psi0[1]) 
  )
)

get_gamma <- function(theta) {
  # theta = -log( (gamma - .3) / (.7 - .3) )
  (parambounds$gamma[2] - parambounds$gamma[1]) * 
    exp(-exp(theta)) + parambounds$gamma[1]
}
get_theta2 <- function(gamma) log(
  -log( 
    (gamma - parambounds$gamma[1]) / (parambounds$gamma[2] - parambounds$gamma[1]) 
  )
)

get_alpha <- function(theta) {
  # theta = log(alpha - 3)
  exp(theta) + parambounds$alpha[1]
}
get_theta3 <- function(alpha) log(alpha - parambounds$alpha[1])

get_beta <- function(theta) {
  # theta = -log( (beta - (-.5)) / (1 - (-.5)) )
  (parambounds$beta[2] - parambounds$beta[1]) * 
    exp(-exp(theta)) + parambounds$beta[1]
}
get_theta4 <- function(beta) log(
  -log( 
    (beta - parambounds$beta[1]) / (parambounds$beta[2] - parambounds$beta[1]) 
  )
)

## Optimization using IPOPT ##
ipopt_objective <- function(theta) ff$fn(theta)
ipopt_objective_gradient <- function(theta) ff$gr(theta)
ipopt_objective_hessian <- function(theta) {
  H <- forceSymmetric(ff$he(theta))
  H <- as(H,"dsTMatrix")
  H
}
ipopt_objective_hessian_x <- function(theta,obj_factor,hessian_lambda) 
  obj_factor * ipopt_objective_hessian(theta)@x
ipopt_objective_hessian_structure <- function(theta) {
  H <- ipopt_objective_hessian(theta)
  H <- as(forceSymmetric(H),'dsTMatrix')
  forStruct = cbind(H@i+1, H@j+1)
  tapply(forStruct[,1], forStruct[,2], c)
}


# Box constraints, to improve stability of optimization
lowerbounds <- c(
  get_theta1(parambounds$Psi0[2] - .001),
  get_theta2(parambounds$gamma[2] - .001),
  get_theta3(parambounds$alpha[1] + .001),
  get_theta4(parambounds$beta[2] - .001)
)

upperbounds <- c(
  get_theta1(parambounds$Psi0[1] + 1),
  get_theta2(parambounds$gamma[1] + .01),
  get_theta3(parambounds$alpha[2] - .01),
  get_theta4(parambounds$beta[1] + .01)
)

# Start in the middle
thetastart <- (upperbounds + lowerbounds)/2


# Nonlinear constraints, specified as a function
ipopt_nonlinear_constraints <- function(theta) Es$fn(theta)

ipopt_nonlinear_constraints_jacobian <- function(theta) {
  J <- Es$gr(theta)
  as(J,"dgTMatrix")
}
ipopt_nonlinear_constraints_jacobian_x <- function(theta) 
  ipopt_nonlinear_constraints_jacobian(theta)@x
ipopt_nonlinear_constraints_jacobian_structure <- function(theta) {
  J <- ipopt_nonlinear_constraints_jacobian(theta)
  J <- as(J,'dgTMatrix')
  forStruct = cbind(J@i+1, J@j+1)
  tapply(forStruct[,2], forStruct[,1], c)
}

nonlinear_lowerbound <- rep(0,nrow(gcdata)+2)
nonlinear_upperbound <- rep(Inf,nrow(gcdata)+2)

# Function and its derivatives
ff <- MakeADFun(data = gcdatalist,
                parameters = list(theta1 = thetastart[1],
                                  theta2 = thetastart[2],
                                  theta3 = thetastart[3],
                                  theta4 = thetastart[4]
                ),
                DLL = "01_astro",
                ADreport = FALSE,
                silent = TRUE)
# Nonlinear constraints and their jacobian
Es <- MakeADFun(data = gcdatalist,
                parameters = list(theta1 = thetastart[1],
                                  theta2 = thetastart[2],
                                  theta3 = thetastart[3],
                                  theta4 = thetastart[4]
                ),
                DLL = "01_astro",
                ADreport = TRUE,
                silent = TRUE)

stopifnot(all(ipopt_nonlinear_constraints(thetastart) > 0))

tm <- Sys.time()
tmp <- capture.output(# Run quietly
  ipopt_result <- ipoptr::ipoptr(
    x0 = (upperbounds + lowerbounds)/2,
    eval_f = ipopt_objective,
    eval_grad_f = ipopt_objective_gradient,
    eval_h = ipopt_objective_hessian_x,
    eval_h_structure = ipopt_objective_hessian_structure(thetastart),
    eval_g = ipopt_nonlinear_constraints,
    eval_jac_g = ipopt_nonlinear_constraints_jacobian_x,
    eval_jac_g_structure = ipopt_nonlinear_constraints_jacobian_structure(thetastart),
    lb = lowerbounds,
    ub = upperbounds,
    constraint_lb = nonlinear_lowerbound,
    constraint_ub = nonlinear_upperbound,
    opts = list(obj_scaling_factor = 1,
                tol = 1e-03)
  )
)
optruntime <- difftime(Sys.time(),tm,units = 'secs')
cat('Run time for mass model optimization:',optruntime,'seconds.\n')

## AGHQ ----
# Create the optimization template
ffa <- list(
  fn = function(x) -1*ff$fn(x),
  gr = function(x) -1*ff$gr(x),
  he = function(x) -1*ff$he(x)
)
useropt <- list(
  ff = ffa,
  mode = ipopt_result$solution,
  hessian = ff$he(ipopt_result$solution)
)
cntrl <- aghq::default_control(negate = TRUE)

# Do the quadrature
tm <- Sys.time()
astroquad <- aghq::aghq(ff,5,thetastart,optresults = useropt,control = cntrl)
quadruntime <- difftime(Sys.time(),tm,units = 'secs')
cat("Run time for mass model quadrature:",quadruntime,"seconds.\n")

# Total time
optruntime + quadruntime # Time difference of 1.010533 secs 2021/04/19

## MCMC ----
stanmod <- tmbstan(
  ff,
  chains = 4,
  cores = 4,
  iter = 1e04,
  warmup = 1e03,
  init = ipopt_result$x0,
  seed = 48645,
  algorithm = "NUTS"
)
# Time
get_elapsed_time(stanmod)
# warmup  sample
# chain:1 1.330630 5.05151
# chain:2 0.912805 4.52549
# chain:3 0.812718 4.90522
# chain:4 0.810394 4.17531
mean(apply(get_elapsed_time(stanmod),1,sum)) # 5.63 seconds

# Inference

# Manually plot
parambounds <- list(
  Psi0 = c(1,200),
  gamma = c(.3,.7), 
  alpha = c(3.0,3.7),
  beta = c(-.5,1)
)

get_psi0 <- function(theta)
  (parambounds$Psi0[2] - parambounds$Psi0[1]) * 
  exp(-exp(theta)) + parambounds$Psi0[1]
get_theta1 <- function(Psi0) 
  log(-log( (Psi0 - parambounds$Psi0[1]) / 
              (parambounds$Psi0[2] - parambounds$Psi0[1]) ))

get_gamma <- function(theta)  
  (parambounds$gamma[2] - parambounds$gamma[1]) * 
  exp(-exp(theta)) + parambounds$gamma[1]
# Add a little buffer, for stability
get_theta2 <- function(gamma) 
  log(-log( (gamma - parambounds$gamma[1] + 1e-03) / 
              (parambounds$gamma[2] - parambounds$gamma[1] + 1e-03) ))

get_alpha <- function(theta)
  exp(theta) + parambounds$alpha[1]
# Add a little buffer, for stability
get_theta3 <- function(alpha) 
  log(alpha - parambounds$alpha[1] + 1e-03)

get_beta <- function(theta)
  (parambounds$beta[2] - parambounds$beta[1]) * 
  exp(-exp(theta)) + parambounds$beta[1]
get_theta4 <- function(beta) 
  log(-log( (beta - parambounds$beta[1]) / 
              (parambounds$beta[2] - parambounds$beta[1]) ))
## Compute the transformed pdfs ##
translist1 <- list(totheta = get_theta1,fromtheta = get_psi0)
translist2 <- list(totheta = get_theta2,fromtheta = get_gamma)
translist3 <- list(totheta = get_theta3,fromtheta = get_alpha)
translist4 <- list(totheta = get_theta4,fromtheta = get_beta)

psi0pdf <- compute_pdf_and_cdf(astroquad$marginals[[1]],translist1)
gammapdf <- compute_pdf_and_cdf(astroquad$marginals[[2]],translist2)
alphapdf <- compute_pdf_and_cdf(astroquad$marginals[[3]],translist3)
betapdf <- compute_pdf_and_cdf(astroquad$marginals[[4]],translist4)

Psi0prior <- function(Psi0) dunif(Psi0,parambounds$Psi0[1],parambounds$Psi0[2],log = FALSE)
gammaprior <- function(gamma) dunif(gamma,parambounds$gamma[1],parambounds$gamma[2],log = FALSE)
alphaprior <- function(alpha) dgamma(alpha - parambounds$alpha[1],shape = 1,rate = 4.6,log = FALSE)
betaprior <- function(beta) dunif(beta,parambounds$beta[1],parambounds$beta[2],log = FALSE)

standata <- as.data.frame(stanmod)
standata$psi0 <- get_psi0(standata[ ,1])
standata$gamma <- get_gamma(standata[ ,2])
standata$alpha <- get_alpha(standata[ ,3])
standata$beta <- get_beta(standata[ ,4])


## Plots ----
PLOTTEXTSIZE <- 28

psi0_postplot <- psi0pdf %>%
  ggplot(aes(x = transparam)) +
  theme_classic() +
  geom_histogram(
    data = standata,
    mapping = aes(x = psi0,y = ..density..),
    bins = 50,
    colour = 'transparent',
    fill = "grey"
  ) +
  geom_line(aes(y = pdf_transparam),size = 1) +
  stat_function(fun = Psi0prior,linetype = "dashed") +
  labs(title = "",x = "",y = "Density") +
  scale_x_continuous(breaks = seq(10,60,by = 5)) +
  coord_cartesian(xlim = c(24,43)) +
  theme(text = element_text(size = PLOTTEXTSIZE))

gamma_postplot <- gammapdf %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_classic() +
  geom_histogram(
    data = standata,
    mapping = aes(x = gamma,y = ..density..),
    bins = 100,
    colour = 'transparent',
    fill = "grey"
  ) +
  geom_line(size = 1) +
  stat_function(fun = gammaprior,linetype = "dashed") +
  labs(title = "",x = "",y = "Density") +
  scale_x_continuous(breaks = seq(0,.44,by = .02)) +
  theme(text = element_text(size = PLOTTEXTSIZE)) +
  coord_cartesian(xlim = c(.3,.4))

alpha_postplot <- alphapdf %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_classic() +
  geom_histogram(
    data = standata,
    mapping = aes(x = alpha,y = ..density..),
    bins = 2000,
    colour = 'transparent',
    fill = "grey"
  ) +
  geom_line(size = 1)  +
  stat_function(fun = alphaprior,linetype = "dashed") +
  labs(title = "",x = "",y = "Density") +
  scale_x_continuous(breaks = seq(3,4,by = .02)) +
  coord_cartesian(xlim = c(3,3.05)) +
  theme(text = element_text(size = PLOTTEXTSIZE))


beta_postplot <- betapdf %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_classic() +
  geom_histogram(
    data = standata,
    mapping = aes(x = beta,y = ..density..),
    bins = 50,
    colour = 'transparent',
    fill = "grey"
  ) +
  geom_line(size = 1)  +
  labs(title = "",x = "",y = "Density") +
  stat_function(fun = betaprior,linetype = "dashed") +
  scale_x_continuous(breaks = seq(-.5,1,by = .1)) +
  coord_cartesian(xlim = c(-.3,.4)) +
  theme(text = element_text(size = PLOTTEXTSIZE))

# Cumulative mass profile
Mr <- function(r,theta) {
  p = get_psi0(theta[1])
  g = get_gamma(theta[2])
  
  # Manual unit conversion into "mass of one trillion suns" (so awesome)
  g*p*r^(1-g) * 2.325e09 * 1e-12
}

rtodo <- 1:150
Mrout <- numeric(length(rtodo))
Mrsdout <- numeric(length(rtodo))
for (rr in 1:length(rtodo)) {
  r <- rtodo[rr]
  
  Mrout[rr] <- compute_moment(
    astroquad$normalized_posterior,
    function(x) Mr(r,x)
  )
  Mrsdout[rr] <- sqrt(compute_moment(
    astroquad$normalized_posterior,
    function(x) (Mr(r,x) - Mrout[rr])^2
  ))
}
cummassplot <- tibble(
  r = rtodo,
  Mr = Mrout,
  Mrsd = Mrsdout
) %>%
  ggplot(aes(x = r)) +
  theme_classic() +
  geom_line(aes(y = Mr),size = 1) +
  geom_ribbon(aes(ymin = Mrout - Mrsdout,
                  ymax = Mrout + Mrsdout),
              alpha = .5) +
  geom_ribbon(aes(ymin = Mrout - 2*Mrsdout,
                  ymax = Mrout + 2*Mrsdout),
              alpha = .2) +
  labs(title = "",x = "",y = bquote('M(r) ('~10^12~M[sun]~')')) +
  scale_x_continuous(breaks = seq(0,150,by = 25)) +
  scale_y_continuous(breaks = seq(0,1,by=.1)) +
  theme(text = element_text(size = PLOTTEXTSIZE))

ggsave(file.path(plotpath,paste0("psi0postplot-",savestamp,".pdf")),psi0_postplot,width = 7,height = 7)
ggsave(file.path(plotpath,paste0("gammapostplot-",savestamp,".pdf")),gamma_postplot,width = 7,height = 7)
ggsave(file.path(plotpath,paste0("alphapostplot-",savestamp,".pdf")),alpha_postplot,width = 7,height = 7)
ggsave(file.path(plotpath,paste0("betapostplot-",savestamp,".pdf")),beta_postplot,width = 7,height = 7)
ggsave(file.path(plotpath,paste0("massplot-",savestamp,".pdf")),cummassplot,width = 7,height = 7)

# KS statistic

aghqpostsamps <- sample_marginal(astroquad,nrow(standata))

suppressWarnings({
  kstable <- data.frame(
    psi0 = ks.test(standata[[1]],aghqpostsamps[[1]])$statistic,
    gamma = ks.test(standata[[2]],aghqpostsamps[[2]])$statistic,
    alpha = ks.test(standata[[3]],aghqpostsamps[[3]])$statistic,
    beta = ks.test(standata[[4]],aghqpostsamps[[4]])$statistic
  )
})

readr::write_csv(kstable,file.path(plotpath,paste0("kstable-",savestamp,".csv")))


