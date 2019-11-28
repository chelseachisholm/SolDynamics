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

coords<- covs[1:2741, 3:4]
distmat <- as.matrix(dist(coords, diag=T, upper=T))
dist_all = abs( distmat - t(distmat))

# Standardize covariates
# mean.elev <- mean(elev, na.rm = TRUE)
# sd.elev <- sd(elev[!is.na(elev)])
# elev <- (elev-mean.elev) / sd.elev
# elev[is.na(elev)] <- 0

## Germ data for sensitivity--< BUT NO HERB DATA< FIX THIS!
germ <- read.table('./data/Germdatanalysis.txt', head=T)

## Herb Chronology data
chrono <- read.csv('./data/Means_of_age.csv', head=T)
chrono <- chrono[1:48, 1:5]
dat <- data.frame(patch=chrono$Name, Elevation=chrono$Elevation, age=chrono$Oldest, y2008=NA, y2009=NA, y2010=NA, y2011=NA, y2012=1) 
dat <- dat %>% mutate(patch=toupper(patch))
dat[which(dat$age == 4), 4:7] <- 1
dat[which(dat$age == 3), 5:7] <- 1
dat[which(dat$age == 2), 6:7] <- 1

#Now match patches to original row locations and add in to z. 
cov$patch <- as.character(cov$Name)
cov$Elevation <- cov$Altitude
rar <-cov %>% reshape2::dcast(PATCH + patch + Elevation~YEAR, value.var='occurrence') %>%
  mutate(patch=toupper(patch))#Issue is named plots are named in patch, but then there is an 'absent' category.
herb<- full_join(rar, dat, by=c('Elevation')) #missing 20 odd patches, not name so have to do it directly... ugh
herb<-herb[,10:14]
#check it worked
colSums(dat[,4:8], na.rm=T)
colSums(herb, na.rm=T) ###YAS QUEEN IT'S THE SAME. No it ain't. For now cut off extras and make sure there are no 1s where there shouldn't be...
herb<- herb[1:2741,]
herb[which(rowSums(occ, na.rm=T)==0),]<-NA
herb[is.na(herb)]<-0

z<-occ

fu <- z+herb
fu[fu==0]<-NA
fu[fu==2]<-1
z<-fu

#Collate data
jags_data <- list(nsite = nrow(occ[1:100,]), nyear = ncol(occ[1:100,]), occ = occ[1:100,], dem = dem[1:100,], flo = flo[1:100,], elev=elev[1:100,], tann=tann[1:100,], D=dist_all[1:100,])
dem_data <- dem
germ_data <- germ





return(jags_data)
}
