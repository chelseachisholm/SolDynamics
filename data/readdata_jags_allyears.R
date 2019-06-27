#### READ IN DATA ####

Importsoldat <- function() {
## Read data
# Bundle data
dat <- read.table('./data/Solidago_PA_2008-2012.txt')
dat %>% group_by(YEAR) %>% filter(occurrence>0) %>% summarize(dem=n())

sols <- dat %>% reshape2::dcast(PATCH~YEAR, value.var='occurrence')
y2009<-c(rep(NA, 2741))
sols$y2009 <- y2009

#need counts in a year x site matrix
y <- as.matrix(sols[,c(2,6,3:5)])
y <- round(y)
y[is.na(y)] <- 0

#read in env data
covs <- read.table('./data/Solidago_PA_2008-2012.txt')

#filter covs to dist_anthro and temperature
elev <- covs  %>% dplyr::select(PATCH, YEAR, Altitude) %>%
  reshape2::dcast(PATCH ~ YEAR, value=Altitude)
elev$y2009 <- elev[,2]
elev <- as.matrix(elev[,c(2,6,3:5)]) #because elev will not change

                
tann <- covs  %>% dplyr::select(PATCH, YEAR, taveSummer) %>%
  reshape2::dcast(PATCH ~ YEAR, value=tave_yy)
tann$y2009 <- tann[,2]
tann <- as.matrix(tann[,c(2,6,3:5)])

                
ddist <- covs  %>% dplyr::select(PATCH, YEAR, Geb) %>% #is Geb correct?
  reshape2::dcast(PATCH ~ YEAR, value=Geb)
ddist$y2009 <- c(rep(NA, 2741))
ddist <- as.matrix(ddist[,c(2,6,3:5)])

coords<- sols[2:3]
pdist <- dist(coords)

# Standardize covariates
# mean.elev <- mean(elev, na.rm = TRUE)
# sd.elev <- sd(elev[!is.na(elev)])
# elev <- (elev-mean.elev) / sd.elev
# elev[is.na(elev)] <- 0

jags_data <- list(nsite = ncol(y), nyear = nrow(y), C = y, year=(c(1:3)-20) / 20, elev=elev, tann=tann, ddist=ddist, D=pdist)

return(jags_data)
}
