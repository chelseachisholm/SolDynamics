### Process demographic data across years 2008-2012 ####
#C. Chisholm, 21.10.2019

library(dplyr)
library(tidyverse)
library(geosphere)
library(rgdal)

sol <- read.table('./data/Solidago_PA_2008-2012.txt')
sol <- sol %>% select(PATCH, YEAR, Easting, Northing, occurrence, Altitude) %>%
  rename(patch=PATCH, year=YEAR, occ=occurrence)  #%>% = éviter l’imbrication des fonctions les unes dans les autres. jardin%>% head() est l'équivalent de head(jardin). Peut être lu comme "ensuite"

sols <- sol %>% spread(year, occ, sep = '') #%>% select(-year2011)
head(sols) #2741 patches x 4 times points (2008,2010,2011,2012)


#### Read in demography data ####
dem8 <- read.table('~/Dropbox/Projects/MIREN/data/To_Gregor_Graubünden_demographic_sampling/insitu_08.txt', header=T )
dem8 <- dem8 %>% mutate(Area=width1*length1, MeanDens = rowMeans(select(., starts_with("d")), na.rm = TRUE)) %>%
  filter(Species=='SC') %>% 
  dplyr::select(Name, Easting, Northing, Area, MeanDens)%>%
  dplyr::rename(Name8=Name)

dem10 <- read.table('~/Dropbox/Projects/MIREN/data/To_Gregor_Graubünden_demographic_sampling/data_in_situ_vor_hint_rhein2010.txt', header=T)
dem10 <- dem10 %>% dplyr::select(-disturbance, -dist_anthro, -DateTime) %>% 
  mutate(Area=width1*length1, MeanDens = rowMeans(select(., paste0("d", 1:20)), na.rm=T),
         MeanHeight= rowMeans(select(., starts_with("Height")), na.rm=T),
         MeanFlower = rowMeans(select(., starts_with("Flower_"))), na.rm=T) %>%
  filter(Species=='SC') %>% 
  dplyr::select(Name, Easting, Northing, Area, MeanDens, MeanFlower, MeanHeight) %>%
  dplyr::rename(Name10=Name)

dem12 <- read.csv('./data/SOLIDAGO_DEMOGRAPHYDATA_50m.csv', sep='\t') #all solidago data
dem12 <- dem12 %>% mutate(MeanHeight = as.numeric(as.character(MeanHeight)), MeanDens=MeanDens, MeanFlower=MeanFlower) %>%
  select(Name, Easting, Northing, AREA, MeanDens, MeanFlower, MeanHeight) %>% dplyr::rename(Area=AREA) %>%
  dplyr::rename(Name12=Name)
#dem12<-read.csv('dem12_updated.csv')

#### Transform UTMs to lat longs ####

#### Join data using space-nearest patch approach ####
#dem8 <- dem8 %>% arrange(Name) %>% mutate(Name=toupper(Name)) #total 122
#dem10 <- dem10 %>% arrange(Name) %>% mutate(Name=toupper(Name)) %>% separate(Name, c('P1','P2'), '-') #2 plots merged, plus a few added, total 131
#dem12 <- dem12 %>% arrange(Name) %>% mutate(Name=toupper(Name)) %>% separate(Name, c('P1','P2','P3','P4'), '-')

ConvertCoordinates <- function(easting,northing) {
  out = cbind(easting,northing)
  mask = !is.na(easting)
  sp <-  sp::spTransform(sp::SpatialPoints(list(easting[mask],northing[mask]),proj4string=sp::CRS("+init=epsg:2056")),sp::CRS("+init=epsg:4326"))
  out[mask,]=sp@coords
  out
}
dem8$latlon <- ConvertCoordinates(dem8$Easting, dem8$Northing)
dem10$latlon <- ConvertCoordinates(dem10$Easting, dem10$Northing)
dem12$latlon <- ConvertCoordinates(dem12$Easting, dem12$Northing)

dem8_nest <- dem8 %>% select(Name8, latlon) #%>% nest(-Name8, Easting, Northing, .key = 'dem8_coords')
dem10_nest <- dem10 %>% select(Name10,  latlon) #%>% nest(-Name10, Easting, Northing, .key = 'dem10_coords')
dem12_nest <- dem12 %>% select(Name12,  latlon) %>% #nest(-Name12, Easting, Northing, .key = 'dem12_coords')



dem.master <- crossing(dem8_nest, dem10_nest)

data.dist <- dem.master %>% 
  mutate(dist = map2_dbl(dem8_coords, dem10_coords, distm)) %>% 
  group_by(Name8) %>% 
  filter(dist == min(dist))

alldem <- full_join(dem8, dem10, by=c('Easting','Northing'))
alldem <- full_join(alldem, dem12, by=c('Easting','Northing'))
alldem <- alldem %>% mutate(PatchArea=sum(AREA.x, AREA.y, AREA, na.rm=T)/3) %>% 
  select(Easting, Northing, PatchArea, MeanDens.x, MeanDens.y, MeanDens, MeanFlower.x, MeanFlower.y, MeanHeight.x, MeanHeight.y) %>%
  rename(Density08 = MeanDens.x, Density10=MeanDens.y, Density12=MeanDens, Flower10 = MeanFlower.x, Flower12=MeanFlower.y, Height10=MeanHeight.x, Height12=MeanHeight.y)

plot(alldem$Density08, alldem$Density10)

plot(alldem$Density10, alldem$Density12)

### JOIN DEMOGRAPHY DATA AND OCCURRENCE DATA


house_nest <- nest(sols, -patch, Easting, Northing, .key = 'patch_coords')
station_nest <- nest(dem, -station_number, .key = 'dem_coords')

data.master <- crossing(house_nest, station_nest)

data.dist <- data.master %>% 
  mutate(dist = map2_dbl(house_coords, station_coords, distm)) %>% 
  group_by(house_number) %>% 
  filter(dist == min(dist))

bla <- full_join(sols, alldem, by=c('Easting', 'Northing'))
rar <- bla %>% arrange(desc(patch)) 
rar <- bla %>% filter(!is.na(patch))
rar$Density10[is.nan(rar$Density10)] <- NA
rar$Density08[is.nan(rar$Density08)] <- NA
#rar$Dens3[is.nan(rar$Dens3)] <- NA

#Create meaningful NA data
#if present and no demo data, NA
#if not present and no demo data, 0 <- FIX THESE
rar$year2008 <- ifelse(!is.na(rar$Density08)&rar$year2008==0, 1, rar$year2008) #no change, 92
rar$year2010 <- ifelse(!is.na(rar$Density10)&rar$year2010==0, 1, rar$year2010)#increased by 2 to 117
rar$year2012 <- ifelse(!is.na(rar$Density12)&rar$year2012==0, 1, rar$year2012) #no chnage, 218

alldem <- rar
write.csv(alldem, './data/cleandem.csv')
# ############
# #Run model of population dynamics
# mt <- glm(occ ~ year, family=binomial, data=sol)
# summary(mt)
# 
# library(car)
# pp <- sum(resid(mt, type='pearson')^2)
# 1-pchisq(pp, mt$df.residual)
# 
# 1-pchisq(mt$deviance, mt$df.residual)
# 
# pp/mt$df.residual
# 
# cr.plots(mt, ask=F)
# influence.measures(mt)
# anova(mt, test="Chisq")
# 
# #plot model of predicted fits
# xs <- seq(2008,2012, l =1000)
# mt.predict  <- predict(mt, type='response', se=T, newdata= data.frame(year = xs))
# plot(occ ~ year, data =sol, xlab='', ylab='', pch=16)
# points(mt.predict$fit ~ xs, type='l', col='grey')
# lines(mt.predict$fit + mt.predict$se.fit ~ xs, col='gray', type='l', lty=2)
# lines(mt.predict$fit - mt.predict$se.fit ~ xs, col='gray', type='l', lty=2)
# mtext('Year', 1, line=3)
# mtext(expression(paste(italic(Solidago), " presence/absence")), 2, line=3)