## 13.4. Analysis of real data set: Single-season occupancy model

library(dplyr)
library(tidyr)
setwd('~/Dropbox/Projects/MIREN/data/For Chelsea 2017')

sol <- read.table('Solidago_PA_2008-2012.txt')
sol <- sol %>% select(PATCH, YEAR, Easting, Northing, occurrence, Altitude) %>%
  rename(patch=PATCH, year=YEAR, occ=occurrence) 
sols <- sol %>% spread(year, occ, sep = '') %>% select(-year2011)
head(sols) #2741 patches x 3 times points (2008,2010,2012)

setwd('~/Desktop/Hohenheim_2017/models/')
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

# Collect the data into suitable structures
y <- as.matrix(sols[, 5:7])
#elev <- sol$Altitude

# Standardize covariates
#mean.elev <- mean(elev, na.rm = TRUE)
#sd.elev <- sd(elev[!is.na(elev)])
#ELEV <- (elev-mean.elev) / sd.elev
#ELEV[is.na(ELEV)] <- 0


#last <- sapply(1:dim(y)[1],
#               function(i) max(grep(FALSE, is.na(y[i, ]))))
y[is.na(y)] <- 0
colnames(y) = NULL
rownames(y) = NULL
stan_data <- list(y = y, nsite = nrow(y), nyear = ncol(y))

## Parameters monitored
params <- c("psi", "phi", "gamma")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
  list(psi1 = runif(1, 0, 1)))

## Call Stan from R
out2 <- stan("solocc.stan",
             data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             control = list(adapt_delta = 0.8),
             open_progress = FALSE)

print(out2, digits = 3)

plot(out2, plotfun = "trace", pars = c("psi"), inc_warmup = TRUE) #all params mixing well
plot(out2, plotfun = "trace", pars = c("gamma"), inc_warmup = TRUE) #all params mixing well
plot(out2, plotfun = "trace", pars = c("phi"), inc_warmup = TRUE) #all params mixing well

plot(out2, pars = c("psi", "phi", "gamma"))

#### INCLUDE DISTANCE ####
coords<- sols[2:3]
D <- as.matrix(dist(coords))
y<-t(y)
stan_data <- list(y = y, nsite = ncol(y), nyear = nrow(y), D=D)

## Parameters monitored
params <- c("psi", "phi", "gamma", "a", "y2")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 1

## Initial values
inits <- lapply(1:nc, function(i)
  list(psi1 = runif(1, 0, 1)))

## Call Stan from R
out2 <- stan("./solocc_distance.stan",
             data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             control = list(adapt_delta = 0.8),
             open_progress = FALSE)

