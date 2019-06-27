library(dplyr)
library(tidyr)

## Read data
# Bundle data
#dat <- read.csv('./data/cleanoccdem.csv')
dat <- read.table('./data/Solidago_PA_2008-2012.txt')
dat %>% group_by(YEAR) %>% filter(occurrence>0) %>% summarize(dem=n())

sols <- dat %>% reshape2::dcast(PATCH~YEAR, value.var='occurrence')
sols$y2009 <- c(rep(NA, 2741))

#need counts in a year x site matrix
y <- as.matrix(sols[,c(2,6,3:5)])
y <- round(y)
y[is.na(y)] <- 0

#read in env data
covs <- read.table('./data/Solidago_PA_2008-2012.txt')

#filter covs to dist_anthro and temperature
elev <- covs  %>% dplyr::select(PATCH, YEAR, Altitude) %>%
  reshape2::dcast(PATCH ~ YEAR, value=Altitude)
elev <- as.matrix(cov1[,c(2:3,4)])

tann <- covs  %>% dplyr::select(PATCH, YEAR, tave_yy) %>%
  reshape2::dcast(PATCH ~ YEAR, value=tave_yy)
tann <- as.matrix(cov2[,c(2:3,4)])

ddist <- covs  %>% dplyr::select(PATCH, YEAR, Geb) %>% #is Geb correct?
  reshape2::dcast(PATCH ~ YEAR, value=Geb)
ddist <- as.matrix(cov2[,c(2:3,4)])

coords<- sols[2:3]
pdist <- dist(coords)

# Standardize covariates
# mean.elev <- mean(elev, na.rm = TRUE)
# sd.elev <- sd(elev[!is.na(elev)])
# elev <- (elev-mean.elev) / sd.elev
# elev[is.na(elev)] <- 0

jags_data <- list(nsite = ncol(y), nyear = nrow(y), C = y, year=(c(1:3)-20) / 20, cov1=cov1, cov2=cov2, D=pdist)

save(jags_data, file="./data/jags_data.RDS")
