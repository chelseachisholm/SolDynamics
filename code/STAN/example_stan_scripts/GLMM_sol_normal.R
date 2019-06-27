## 4. Poisson GLMM for random effects
## 4.2. Accounting for overdispersion by random effects-modeling in R
## and WinBUGS
## 4.2.1. Generation and analysis of simulated data
setwd('~/Desktop/Hohenheim_2017/data/')
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Read data
# Bundle data
dat <- read.csv('cleanoccdem.csv')
sols <- dat %>% gather(patch) %>% spread(year, Dens)
sols <- dat %>% dcast(year~patch, value.var='Dens', fun.aggregate = sum)

#need counts in a year x site matrix
y <- as.matrix(sols[, 2:2742])
y <- round(y)
y[is.na(y)] <- 0

stan_data <- list(nsite = ncol(y), nyear = nrow(y), C = y, year=(c(1:3)-20) / 20)

# Initial values
inits <- function() list(mu = runif(1, 0, 2), alpha = runif(stan_data$nsite, -1, 1), beta = runif(2, -1, 1), sd.alpha = runif(1, 0, 0.1), sd.year = runif(1, 0, 0.1))

# Parameters monitored (may want to add "lambda")
params <- c("mu", "alpha", "beta", "sd_alpha", "sd_year")

# MCMC settings (may have to adapt)
ni <- 10000
nt <- 10
nb <- 5000
nc <- 4


## Call Stan from R
out <- stan("../GLMM_sol_normal_linear.stan", data = stan_data,
            init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)

#Stan does not support NA data? Maybe if I add 0s?
  ## Summarize posteriors
print(out)
pairs(out)
plot(out, plotfun = "trace", pars = c("beta"), inc_warmup = TRUE) #meh, but also one chain
plot(out, pars = c("beta"))
