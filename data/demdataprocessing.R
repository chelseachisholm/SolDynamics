### Process demographic data across years 2008-2012 ####
#C. Chisholm, 21.10.2019

library(dplyr)
library(tidyr)
setwd('./')

sol <- read.table('Solidago_PA_2008-2012.txt')
sol <- sol %>% select(PATCH, YEAR, Easting, Northing, occurrence, Altitude) %>%
  rename(patch=PATCH, year=YEAR, occ=occurrence)  #%>% = éviter l’imbrication des fonctions les unes dans les autres. jardin%>% head() est l'équivalent de head(jardin). Peut être lu comme "ensuite"

sols <- sol %>% spread(year, occ, sep = '') #%>% select(-year2011)
head(sols) #2741 patches x 4 times points (2008,2010,2011,2012)


############

#Read in density dependence
dem8 <- read.table('~/Dropbox/Projects/MIREN/data/To_Gregor_Graubünden_demographic_sampling/insitu_08.txt', header=T )
dem8 <- dem8 %>% mutate(AREA=width1*length1, MeanDens = rowMeans(select(., starts_with("d")), na.rm = TRUE)) %>%
  filter(Species=='SC') %>% 
  select(Name, Easting, Northing, AREA, MeanDens)

dem10 <- read.table('~/Dropbox/Projects/MIREN/data/To_Gregor_Graubünden_demographic_sampling/data_in_situ_vor_hint_rhein2010.txt', header=T)
dem10 <- dem10 %>% select(-c(disturbance, dist_anthro, DateTime)) %>% 
  mutate(AREA=width1*length1, MeanDens = rowMeans(select(., paste0("d", 1:20)), na.rm=T),
         MeanHeight= rowMeans(select(., starts_with("Height")), na.rm=T),
         MeanFlower = rowMeans(select(., starts_with("Flower_"))), na.rm=T) %>%
  filter(Species=='SC') %>% 
  select(Name, Easting, Northing, AREA, MeanDens, MeanFlower, MeanHeight)

dem12 <- read.csv('SOLIDAGO_DEMOGRAPHYDATA_50m.csv', sep='\t') #all solidago data
dem12 <- dem12 %>% mutate(MeanHeight = as.numeric(as.character(MeanHeight)), MeanDens=MeanDens, MeanFlower=MeanFlower) %>%
  select(Name, Easting, Northing, AREA, MeanDens, MeanFlower, MeanHeight)
#dem12<-read.csv('dem12_updated.csv')

############
#Assess density dependence
#Using full join
dem8 <- dem8 %>% arrange(Name) %>% mutate(Name=toupper(Name)) #total 122
dem10 <- dem10 %>% arrange(Name) %>% mutate(Name=toupper(Name)) %>% separate(Name, c('P1','P2'), '-') #2 plots merged, plus a few added, total 131
dem12 <- dem12 %>% arrange(Name) %>% mutate(Name=toupper(Name)) %>% separate(Name, c('P1','P2','P3','P4'), '-')
alldem <- full_join(dem8, dem10, by=c('Easting','Northing'))
alldem <- full_join(alldem, dem12, by=c('Easting','Northing'))
alldem <- alldem %>% select(Easting, Northing, MeanDens.x, MeanDens.y, MeanDens, MeanFlower.x, MeanFlower.y, MeanHeight.x, MeanHeight.y) %>%
  rename(Density08 = MeanDens.x, Density10=MeanDens.y, Density12=MeanDens, Flower10 = MeanFlower.x, Flower12=MeanFlower.y, Height10=MeanHeight.x, Height12=MeanHeight.y)

plot(alldem$Density08, alldem$Density10)

plot(alldem$Density10, alldem$Density12)

### JOIN DEMOGRAPHY DATA AND OCCURRENCE DATA
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
write.csv(alldem, 'cleandem.csv')
############
#Run model of population dynamics
mt <- glm(occ ~ year, family=binomial, data=sol)
summary(mt)

library(car)
pp <- sum(resid(mt, type='pearson')^2)
1-pchisq(pp, mt$df.residual)

1-pchisq(mt$deviance, mt$df.residual)

pp/mt$df.residual

cr.plots(mt, ask=F)
influence.measures(mt)
anova(mt, test="Chisq")

#plot model of predicted fits
xs <- seq(2008,2012, l =1000)
mt.predict  <- predict(mt, type='response', se=T, newdata= data.frame(year = xs))
plot(occ ~ year, data =sol, xlab='', ylab='', pch=16)
points(mt.predict$fit ~ xs, type='l', col='grey')
lines(mt.predict$fit + mt.predict$se.fit ~ xs, col='gray', type='l', lty=2)
lines(mt.predict$fit - mt.predict$se.fit ~ xs, col='gray', type='l', lty=2)
mtext('Year', 1, line=3)
mtext(expression(paste(italic(Solidago), " presence/absence")), 2, line=3)