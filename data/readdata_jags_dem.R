#### READ IN DATA ####

Importsoldat <- function() {

#### SOLIDAGO DATA ####

# Occurrence data
sol <- read.table('./data/Solidago_PA_2008-2012.txt')
sol <- sol %>% select(PATCH, YEAR, Easting, Northing, occurrence, Altitude) %>%
    rename(patch=PATCH, year=YEAR, occ=occurrence)  #%>% = éviter l’imbrication des fonctions les unes dans les autres. jardin%>% head() est l'équivalent de head(jardin). Peut être lu comme "ensuite"
  
sols <- sol %>% spread(year, occ, sep = '') #%>% select(-year2011)
head(sols) #2741 patches x 4 times points (2008,2010,2011,2012)
cov <- read.table('./data/Solidago_PA_2008-2012.txt')
# sols <-cov %>% reshape2::dcast(PATCH~YEAR, value.var='occurrence') %>% rename(patch=PATCH)

# Demography data
dat <- read.csv('./data/cleandem.csv')
y2009<-c(rep(NA, 2741))
dat$missing <- y2009
occ <- dat[,c('year2008', 'missing', 'year2010', 'year2011', 'year2012')] %>% as.matrix(.)
dem <- dat[,c('Density08', 'missing', 'Density10','missing','Density12')] %>% as.matrix(.)
flo <- dat[,c('missing', 'missing', 'Flower10', 'missing', 'Flower12')] %>% as.matrix(.)
hei <- dat[,c('missing', 'missing', 'Height10', 'missing', 'Height12')] %>% as.matrix(.)
pat <-dat[,c('Area08', 'missing', 'Area10', 'missing', 'Area12')] %>% as.matrix(.)

#### COVARIATES ####
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

gdd <- covs  %>% dplyr::select(PATCH, YEAR, ddeg000_h) %>%
  reshape2::dcast(PATCH ~ YEAR, value=taveSummer)
gdd$y2009 <- gdd[,2]
gdd <- as.matrix(gdd[,c(2:6)])

fdd <- covs  %>% dplyr::select(PATCH, YEAR, sfroyy_100) %>%
  reshape2::dcast(PATCH ~ YEAR, value=taveSummer)
fdd$y2009 <- fdd[,2]
fdd <- as.matrix(fdd[,c(2:6)])

coords<- covs[1:2741, 3:4]
distmat <- as.matrix(dist(coords, diag=T, upper=T))
dist_all = abs( distmat - t(distmat))


#### DETECTION DATA ####

# Herb Chronology data (generate presence/absences)
chrono <- read.csv('./data/Means_of_age_matched.csv', head=T) #48 rows, has matched columns from germ dataset
chrono <- chrono %>% mutate(Altitude=(Elevation +Elevation.1)/2)
dat <- data.frame(Name=chrono$Name, Elevation=chrono$Altitude, Easting=chrono$Easting, Northing=chrono$Northing, age=chrono$Oldest, y2008=NA, y2009=NA, y2010=NA, y2011=NA, y2012=1) 
dat[which(dat$age == 4), 6:9] <- 1
dat[which(dat$age == 3), 7:9] <- 1
dat[which(dat$age == 2), 8:9] <- 1

#Match to original pa data
ConvertCoordinates <- function(easting,northing) {
  out = cbind(easting,northing)
  mask = !is.na(easting)
  sp <-  sp::spTransform(sp::SpatialPoints(list(easting[mask],northing[mask]),proj4string=sp::CRS("+init=epsg:2056")),sp::CRS("+init=epsg:4326"))
  out[mask,]=sp@coords
  out
}

dat$latlon <- ConvertCoordinates(dat$Easting, dat$Northing)

#Convert original datast UTMNs
pa <- sols
pa$latlon <- ConvertCoordinates(pa$Easting, pa$Northing)

dat_nest <- dat %>% select(Name, latlon) %>% nest(-Name, latlon, .key = 'dat_coords')
pa_nest <- pa %>% select (patch, latlon) %>% nest(-patch, latlon, .key = 'pa_coords')

#### Join data using space-nearest patch approach ####
dem.master <- crossing(dat_nest, pa_nest)

herb_d <- dem.master %>% 
  mutate(dist = map2_dbl(dat_coords, pa_coords, distm)) %>% 
  group_by(Name) %>% 
  filter(dist == min(dist))

# Match patches using join
sols <- sols %>% rename(Elevation=Altitude)
herbchrono <- full_join(dat, herb_d, by='Name') %>% 
  filter(dist<50) %>% #loses 12 patches (making this 100 removes 6)
  select(patch, y2008, y2009, y2010, y2011, y2012) 
herbchrono[duplicated(herbchrono$patch),] # some patches are combined
herbchrono <- herbchrono[-c(20:21,31,44),]
herbchrono[is.na(herbchrono)]<-0

herb<- full_join(sols, herbchrono) #%>% arrange(patch)  #missing 8 patches, greater than 200 m away from track grids
#There appears to be 16 patches missing P.1400, P/1408, P.1412, P.1416. Others missing later on as well. Why?

# Check matched correctly

dim(herb[which(rowSums(herb[,9:13])>0),]) #33 patches, that's right! But it looks like there are some false positives in our dataset. Remove for now.
herb[which(rowSums(herb[,9:13])>0),]

# Remove false positives
herbs <- herb %>% mutate(y2008=ifelse(year2008==1 & y2008==0, NA, y2008), 
                y2010=ifelse(year2010==1 & y2010==0, NA, y2010),
                y2011=ifelse(year2011==1 & y2011==0, NA, y2011),
                y2012=ifelse(year2012==1 & y2012==0, NA, y2012))

herb <- herbs[,c(9:13)]
herb[is.na(herb)] <- 0
rar<-occ
rar[is.na(rar)] <- 0

fu <- rar+herb
fu[fu==0]<-NA
fu[fu==2]<-1
z<-fu

#Collate data
#occ[seq(1, nrow(occ), 50),]
jags_data <- list(nsite = nrow(occ), nyear = ncol(occ), occ = occ, z = z, dem = dem, flo = flo, pat = pat, elev=elev, distmat=distmat)

save(jags_data, file='jags_data.RData')

return(jags_data)
}
