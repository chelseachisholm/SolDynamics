#### READ IN DATA ####

Importsoldat <- function() {
## Read data
# occurrence data 
cov <- read.table('./data/Solidago_PA_2008-2012.txt')
sols <-cov %>% reshape2::dcast(PATCH~YEAR, value.var='occurrence') %>% rename(patch=PATCH)

# Demography data
dat <- read.csv('./data/cleandem.csv')
dim(dat) #doesn't match sols up above!
dat %>%
  filter(duplicated(patch) | duplicated(patch, fromLast = TRUE))
dat <-dat[-65,]
y2009<-c(rep(NA, 2741))
dat$y2009 <- y2009
occ <- as.matrix(dat[,c(6, 17, 7, 8, 9)])
dem <- as.matrix(dat[,c(10,17,11,17, 12)])
flo <- as.matrix(dat[,c(17,17,13,17, 14)])
hei <- as.matrix(dat[,c(17,17,15,17, 16)])

## Add in covariates
#read in env data
covs <- read.table('./data/Solidago_PA_2008-2012.txt')

#filter covs to dist_anthro and temperature
elev <- covs  %>% dplyr::select(PATCH, YEAR, Altitude) %>%
  reshape2::dcast(PATCH ~ YEAR, value=Altitude) %>% rename(patch=PATCH)
elev$y2009 <- elev[,2]
elev <- as.matrix(elev[,c(2,6,3:5)]) #because elev will not change

tann <- covs  %>% dplyr::select(PATCH, YEAR, taveSummer) %>%
  reshape2::dcast(PATCH ~ YEAR, value=taveSummer)
tann$y2009 <- tann[,2]
tann <- as.matrix(tann[,c(2,6,3:5)])

#ddist <- covs  %>% dplyr::select(PATCH, YEAR, Geb) %>% #is Geb correct?
#  reshape2::dcast(PATCH ~ YEAR, value=Geb)
#ddist$y2009 <- c(rep(NA, 2741))
#ddist <- as.matrix(ddist[,c(2,6,3:5)])

coords<- covs[3:4]
distmat <- as.matrix(dist(coords, diag=T, upper=T))
dist_all = abs( distmat - t(distmat))

# Standardize covariates
# mean.elev <- mean(elev, na.rm = TRUE)
# sd.elev <- sd(elev[!is.na(elev)])
# elev <- (elev-mean.elev) / sd.elev
# elev[is.na(elev)] <- 0

jags_data <- list(nsite = nrow(occ), nyear = ncol(occ), occ = occ, dem = dem, flo = flo, hei = hei, elev=elev, tann=tann, D=dist_all)

return(jags_data)
}
