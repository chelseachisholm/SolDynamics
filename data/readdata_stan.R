library(dplyr)
library(tidyr)

## Read data
# Bundle data
dat <- read.csv('./data/cleanoccdem.csv')
dat %>% group_by(year) %>% filter(Dens>0) %>% summarize(dem=n())

sols <- dat %>% reshape2::dcast(year~patch, value.var='occ', fun.aggregate = sum)

#need counts in a year x site matrix
y <- as.matrix(sols[, 2:2742])
y <- round(y)
y[is.na(y)] <- 0

#read in env data
covs <- read.table('./data/Solidago_PA_2008-2012.txt')

#filter covs to dist_anthro and temperature
cov1 <- covs  %>% dplyr::select(PATCH, YEAR, Altitude) %>%
  reshape2::dcast(PATCH ~ YEAR, value=Altitude)
cov1 <- as.matrix(cov1[, 2:2742])
y <- round(y)
y[is.na(y)] <- 0

cov2 <- covs  %>% dplyr::select(PATCH, YEAR, taveSummer) %>%
  reshape2::dcast(PATCH ~ YEAR, value=taveSummer)

stan_data <- list(nsite = ncol(y), nyear = nrow(y), C = y, year=(c(1:3)-20) / 20, cov1)

