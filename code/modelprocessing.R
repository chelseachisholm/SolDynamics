#### Process Output from JAGS models ####
#07.06.2019 C. Chisholm

require(mcmcplots)
library(MCMCvis)
library(ggplot2)
library(cowplot)

source('~/Dropbox/Projects/Invasion Dynamics/models/code/helperfunctions.R')
#load("jags_data.RData")
load('~/Dropbox/projects/Invasion Dynamics/models/modeloutput122020.RData')
Samples <- parsamples$mcmc
mcmc <- as.data.frame(cbind(Samples[[1]], Samples[[2]], Samples[[3]]))
#load('~/Dropbox/projects/Invasion Dynamics/models/jags_data.RData')


### Convergence checks ####
# Plot the mcmc chain and the posterior sample for p
plot(Samples, ask=T)

# convergence check
gelman.plot(Samples)
coda::gelman.diag(Samples)

# Statistical summaries of the posterior sample
summary(Samples)

Samples.wd <- window(Samples, start = 1500)
plot(Samples.wd)

caterplot(Samples, "sigma_flo", style = "plain")

n_occ <- rev(c(98, NA, 123, 135, 224, NA, NA, NA, NA, NA, NA, NA, 620))
caterplot(Samples, "n_occ", style = "plain")
caterpoints(n_occ) #Can add the true occupancy here!
coefs = mcmc[, c("n_occ[1]", "n_occ[2]", "n_occ[3]", "n_occ[4]", "n_occ[5]", "n_occ[13]")]
colMeans(coefs)
sd(coefs)



MCMCsummary(Samples, params = 'beta_phi')
MCMCsummary(Samples, params = 'beta_rho')

#### Plot: p and alpha ####
coefs = mcmc[, "p"]
dat <- reshape2::melt(coefs)

ggplot(dat, aes(x = value)) +
  geom_histogram() +
  theme(legend.position = "none") +
  labs(title = 'Model estimates') +
  scale_y_discrete("Parameter p") +
  scale_x_continuous("Value") + 
  theme_classic(base_size = 16) +
  scale_fill_manual(values="seagreen") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )

coefs = mcmc[, "log.alpha"]
dat <- reshape2::melt(coefs)
MCMCsummary(Samples, params = 'log.alpha')
caterplot(parsamples$mcmc, "log.alpha", style = "plain")

ggplot(dat, aes(x = value)) +
  geom_histogram() +
  theme(legend.position = "none") +
  labs(title = 'Alpha') +
  scale_y_discrete("Parameter alpha") +
  scale_x_continuous("Value") + 
  theme_classic(base_size = 16) +
  scale_fill_manual(values="seagreen") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )



#### Plot 1: Occupied patches ####
coefs = mcmc[, c("n_occ[1]", "n_occ[2]", "n_occ[3]", "n_occ[4]", "n_occ[5]")] #,
#"n_occ[13]")]
dat <- reshape2::melt(coefs)
dat$variable <- factor(dat$variable, labels=c('2008', '2009', '2010', '2011', '2012',
                                              '2020'))

# library
library(ggridges)

# basic example
n_occ <- data.frame(year=c('2008','2009','2010','2011','2012','2020'), nocc=c(98, NA, 123, 135, 224, 620))
ggplot(dat, aes(x = value, y = variable, fill="#69b3a2")) +
  geom_boxplot() +
  #theme_ridges() + 
  theme(legend.position = "none") +
  labs(title = 'Total Patch Occupancy') +
  scale_y_discrete("Year") +
  scale_x_continuous("No. Patches Occupied") + 
  theme_classic(base_size = 16) +
  scale_fill_manual(values="seagreen") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  geom_point(data=n_occ, aes(x=nocc, y=year))

#### Plot 2: Overall C, E and Turnover over the years ####
#Extinction
coefs = mcmc[, c("sext[2]", "sext[3]", "sext[4]", "sext[5]")]
dat <- reshape2::melt(coefs)
dat$variable <- factor(dat$variable, labels=c('2009', '2010', '2011', '2012'))

# basic example
pext <- ggplot(dat, aes(x = value, y = variable, fill="#69b3a2", alpha=0.7)) +
  geom_density_ridges() +
  #theme_ridges() + 
  theme(legend.position = "none") +
  scale_y_discrete("Year") +
  scale_x_continuous("No. Patches") + 
  theme_classic(base_size = 16) +
  scale_fill_manual(values='#f2ad73') +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )
pext

#Colonisation
coefs = mcmc[, c("scol[2]", "scol[3]", "scol[4]", "scol[5]")]
dat <- reshape2::melt(coefs)
dat$variable <- factor(dat$variable, labels=c('2009', '2010', '2011', '2012'))

# basic example
pcol <- ggplot(dat, aes(x = value, y = variable, fill="#69b3a2", alpha=0.7)) +
  geom_density_ridges() +
  #theme_ridges() + 
  theme(legend.position = "none") +
  scale_y_discrete("Year") +
  scale_x_continuous("No. Patches") + 
  theme_classic(base_size = 16) +
  scale_fill_manual(values='#f2ad73') +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )
pcol

#Turnover
coefs = mcmc[, c("stur[2]", "stur[3]", "stur[4]", "stur[5]")]
dat <- reshape2::melt(coefs)
dat$variable <- factor(dat$variable, labels=c('2009', '2010', '2011', '2012'))

# basic example
ptur <- ggplot(dat, aes(x = value, y = variable, fill="#69b3a2", alpha=0.7)) +
  geom_density_ridges() +
  #theme_ridges() + 
  theme(legend.position = "none") +
  scale_y_discrete("Year") +
  scale_x_continuous("No. Patches") + 
  theme_classic(base_size = 16) +
  scale_fill_manual(values='#f2ad73') +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )
ptur

plot_grid(pcol, pext, ptur, nrow=1)

#### Plot 1: Survival  ####
elev <- jags_data$elev
nelev<- scale(elev)

occ <- jags_data$occ[,4]
df <- data.frame(elev=elev, occ=occ)

## Calculate the fitted values
nvalues <- 100
newdata <- seq(min(nelev), max(nelev), length.out = nvalues)

pred_mean_dist <- matrix(NA, nrow = nrow(mcmc), ncol = nvalues)

for (i in 1:nrow(pred_mean_dist)){
  pred_mean_dist[i,] <- antilogit(mcmc[i,"beta_phi[1]"] + newdata * mcmc[i,"beta_phi[2]"])
}

credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)

pred_y <- antilogit(mean(mcmc[,'beta_phi[1]']) + newdata* mean(mcmc[,'beta_phi[2]']))
newdata_x <- seq(min(elev), max(elev), length.out = nvalues)
dat <- data.frame(elev=newdata_x, estimate=pred_y, l_cred=credible_lower, u_cred=credible_upper)

p1 <- ggplot(dat, aes(y = estimate, x = elev)) + geom_line() + #estimate is P(extinction)
  geom_ribbon(aes(ymin = l_cred, ymax = u_cred), fill = "orange", alpha = 0.3) + 
  #geom_count(data = df, mapping = aes(x = elev, y = occ), fill='black', alpha=0.3) +
  theme(legend.position="none")  +
  #scale_size_area(max_size = 10) +
  scale_y_continuous("P(Survival)") +
  scale_x_continuous("Elevation (m)") + theme_classic(base_size = 16)

p1

#### Plot 2: Dispersal  ####

nelev<- scale(elev)

alpha <- mcmc[, grep("alpha", colnames(mcmc), fixed=T)]

## Calculate the fitted values

alpha_mean <- colMeans(alpha)
cl_1 <- apply(alpha, MARGIN = 2, quantile, prob = 0.025)
cu_1 <- apply(alpha, MARGIN = 2, quantile, prob = 0.975)

dat <- data.frame(estimate=alpha_mean, l_cred=cl_1, u_cred=cu_1)

p1 <- ggplot(dat, aes(y = estimate)) + geom_smooth() + 
  geom_ribbon(aes(ymin = (l_cred), ymax = (u_cred)), fill = "orange", alpha = 0.3) + 
  scale_y_continuous("P(Colonisation)") +
  scale_x_continuous("Elevation (m)") + theme_classic(base_size = 16)

p1 #This is a bit fucked- it isn't with respect to elevation, so I'm not sure how to show it in a smooth fashion (vs. the dumb caterplot I have now)
#Maybe gamma0 ~ elevation?

#### Plot 3: Population size ####

# Calculate max elevatio of occupancy (from df in Plot 1)
dfm <- df[df$occ==1,]
dfm[which.max(dfm$elev),] #1112

## Calculate the fitted values
nvalues <- 100
newdata <- seq(min(nelev), max(nelev), length.out = nvalues)

pred_mean_dist <- matrix(NA, nrow = nrow(mcmc), ncol = nvalues)

for (i in 1:nrow(pred_mean_dist)){
  pred_mean_dist[i,] <- (mcmc[i,"beta_N[1]"] + newdata * mcmc[i,"beta_N[2]"])
}

credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)

pred_y <- (mean(mcmc[,'beta_N[1]']) + newdata* mean(mcmc[,'beta_N[2]']))
newdata_x <- seq(min(elev), max(elev), length.out = nvalues)
dat <- data.frame(elev=newdata_x, estimate=pred_y, l_cred=credible_lower, u_cred=credible_upper)

p3 <- ggplot(dat, aes(y = exp(estimate), x = elev)) + geom_line() +
  geom_ribbon(aes(ymin = exp(l_cred), ymax = exp(u_cred)), fill = "orange", alpha = 0.3) +
  scale_y_continuous("Abundance") +
  scale_x_continuous("") + theme_classic(base_size = 16) + 
  geom_vline(xintercept = 1112, linetype="dashed", color = "darkgrey", size=1.5)

p3
# p3 <- ggplot(dat, aes(y = (estimate*15+15), x = elev)) + geom_line() + 
#   geom_ribbon(aes(ymin = (l_cred*15+15), ymax = (u_cred*15+15)), fill = "orange", alpha = 0.3) + 
#   scale_y_continuous("Abundance") +
#   scale_x_continuous("") + theme_classic(base_size = 16)

#### Plot 4: Flowering ####

## Calculate the fitted values
nvalues <- 100
newdata <- seq(min(nelev), max(nelev), length.out = nvalues)

pred_mean_dist <- matrix(NA, nrow = nrow(mcmc), ncol = nvalues)

for (i in 1:nrow(pred_mean_dist)){
  pred_mean_dist[i,] <- (mcmc[i,"beta_flo[1]"] + newdata * mcmc[i,"beta_flo[2]"])
}

credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)

pred_y <- (mean(mcmc[,'beta_flo[1]']) + newdata* mean(mcmc[,'beta_flo[2]']))
newdata_x <- seq(min(elev), max(elev), length.out = nvalues)
dat <- data.frame(elev=newdata_x, estimate=pred_y, l_cred=credible_lower, u_cred=credible_upper)

p4 <- ggplot(dat, aes(y = exp(estimate), x = elev)) + geom_line() + 
  geom_ribbon(aes(ymin = exp(l_cred), ymax = exp(u_cred)), fill = "orange", alpha = 0.3) + 
  scale_y_continuous("Total No. Flowers") +
  scale_x_continuous("Elevation (m)") + theme_classic(base_size = 16) + 
  geom_vline(xintercept = 1112, linetype="dashed", color = "darkgrey", size=1.5)

p4


#### Plot for all demographic parameters ####
library("cowplot")
plot_grid(p3, p4, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

#### Plot 5 & 6: Survival and Colonisation ####
## Calculate the fitted values
nvalues <- 100
newdata <- seq(min(nelev), max(nelev), length.out = nvalues)

pred_mean_dist <- matrix(NA, nrow = nrow(mcmc), ncol = nvalues)

for (i in 1:nrow(pred_mean_dist)){
  pred_mean_dist[i,] <- antilogit(mcmc[i,"beta_phi[1]"] + newdata * mcmc[i,"beta_phi[2]"])
}

credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)

pred_y <- antilogit(mean(mcmc[,'beta_phi[1]']) + newdata* mean(mcmc[,'beta_phi[2]']))
newdata_x <- seq(min(elev), max(elev), length.out = nvalues)
dat5 <- data.frame(elev=newdata_x, estimate=pred_y, l_cred=credible_lower, u_cred=credible_upper)

p5 <- ggplot(dat5, aes(y = estimate, x = elev)) + geom_line() + 
  geom_ribbon(aes(ymin = l_cred, ymax = u_cred), fill = "orange", alpha = 0.3) + 
  scale_y_continuous("P(Survival)") +
  scale_x_continuous("Elevation (m)") + theme_classic(base_size = 16) + 
  geom_vline(xintercept = 1112, linetype="dashed", color = "darkgrey", size=1.5)

p5

## Calculate the fitted values
nvalues <- 100
newdata <- seq(min(nelev), max(nelev), length.out = nvalues)

pred_mean_dist <- matrix(NA, nrow = nrow(mcmc), ncol = nvalues)

for (i in 1:nrow(pred_mean_dist)){
  pred_mean_dist[i,] <- antilogit(mcmc[i,"beta_rho[1]"] + newdata * mcmc[i,"beta_rho[2]"])
}

credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)

pred_y <- antilogit(mean(mcmc[,'beta_rho[1]']) + newdata* mean(mcmc[,'beta_rho[2]']))
newdata_x <- seq(min(elev), max(elev), length.out = nvalues)
dat6 <- data.frame(elev=newdata_x, estimate=pred_y, l_cred=credible_lower, u_cred=credible_upper)

p6 <- ggplot(dat6, aes(y = estimate, x = elev)) + geom_line() + 
  geom_ribbon(aes(ymin = l_cred, ymax = u_cred), fill = "orange", alpha = 0.3) + 
  scale_y_continuous("P(Colonisation)") +
  scale_x_continuous("Elevation (m)") + theme_classic(base_size = 16) + 
  geom_vline(xintercept = 1112, linetype="dashed", color = "darkgrey", size=1.5)

p6


#### Plot for all colonisation parameters ####
library("cowplot")
plot_grid(p5, p6,  
          #labels = c("A", "B"),
          ncol = 2, nrow = 1)

#### Plot 7: Combining colonisation and extinction rates (BUT NEED TO EXTRACT GAMMA!) #####

# ## Calculate the fitted values
# nvalues <- 100
# newdata <- seq(min(nelev), max(nelev), length.out = nvalues)
# 
# newdata2 <- rep(mean(anthro), length=nvalues)
# 
# pred_mean_dist <- matrix(NA, nrow = nrow(mcmc), ncol = nvalues)
# 
# for (i in 1:nrow(pred_mean_dist)){
#   pred_mean_dist[i,] <- antilogit(mcmc[i,"beta_rho[1]"] + newdata * mcmc[i,"beta_rho[2]"] + newdata2 * mcmc[i,"beta_rho[3]"])
# }
# 
# credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
# credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
# 
# pred_y <- antilogit(mean(mcmc[,'beta_rho[1]']) + newdata* mean(mcmc[,'beta_rho[2]']) + newdata2 * mcmc[i,"beta_rho[3]"])
# newdata_x <- seq(min(elev), max(elev), length.out = nvalues)
# dat1 <- data.frame(elev=newdata_x, estimate=pred_y, l_cred=credible_lower, u_cred=credible_upper)
# 
# #### Extinction
# ## Calculate the fitted values
# nvalues <- 100
# newdata <- seq(min(nelev), max(nelev), length.out = nvalues)
# 
# pred_mean_dist <- matrix(NA, nrow = nrow(mcmc), ncol = nvalues)
# 
# for (i in 1:nrow(pred_mean_dist)){
#   pred_mean_dist[i,] <- antilogit(mcmc[i,"beta_phi[1]"] + newdata * mcmc[i,"beta_phi[2]"])
# }
# 
# credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
# credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
# 
# pred_y <- antilogit(mean(mcmc[,'beta_phi[1]']) + newdata* mean(mcmc[,'beta_phi[2]']))
# newdata_x <- seq(min(elev), max(elev), length.out = nvalues)
# dat2 <- data.frame(elev=newdata_x, estimate=pred_y, l_cred=credible_lower, u_cred=credible_upper)
# 
# dat <- bind_rows('Colonisation' = dat1, 'Extinction' = dat2, .id = 'process')
# 
# p7 <- ggplot(dat, aes(y = estimate, x = elev, color=process)) + geom_line() + #estimate is P(extinction)
#   geom_ribbon(aes(ymin = (l_cred), ymax = (u_cred), fill = process), alpha = 0.3) +
#   theme(legend.position="none")  +
#   scale_y_continuous("P(process)") +
#   scale_x_continuous("Elevation (m)") + theme_classic(base_size = 16)
# 
# p7
# 
# p7 <- ggplot(dat, aes(y = estimate, x = elev)) + geom_line() + 
#   geom_ribbon(aes(ymin = l_cred, ymax = u_cred), fill = "orange", alpha = 0.3) + 
#   scale_y_continuous("Colonisation") +
#   scale_x_continuous("Elevation (m)") + theme_classic(base_size = 16) + 
#   geom_vline(xintercept = 1112, linetype="dashed", color = "darkgrey", size=1.5)
# 
# p7
#### EXTRA CODE ####

# Plot landscape connectivity 

par_vals = summary(Samples)$statistics
gammas <- par_vals %>% as.data.frame() %>% rownames_to_column(var='param') %>%
  filter(grepl('gamma', param)) %>% filter(param != 'gamma0') %>%
  mutate(Year= rep(1:4, each=2741), Patch=rep(1:2741, times=4), Elev=rep(elev[,1], times=4))

gammas %>% group_by(Year) %>% summarize(mean=mean(Mean))
ggplot(gammas, aes(x=Elev, y=Mean)) + geom_point() + facet_wrap(~Year) 


gamma0 <- par_vals %>% as.data.frame() %>% rownames_to_column(var='param') %>%
  filter(grepl('gamma', param)) %>% filter(param == 'gamma0') 

#Using simulated data
range(elev)
# range(flo, na.rm=T)
# range(dem, na.rm=T)

nelev<- scale(elev)
mcmc <- as.data.frame(cbind(Samples[[1]], Samples[[2]], Samples[[3]])) 
nvalues <- 100
newdata <- seq(min(nelev), max(nelev), length.out = nvalues)
p_surv = antilogit(par_vals[1,1] + par_vals[2,1]*newdata)
# p_flo = antilogit(par_vals[3,1] + par_vals[4,1]*c(0:46))
# p_abu = antilogit(par_vals[1,1] + par_vals[2,1]*c(0:70))

sim_surv <- data.frame(Elev=sample(elev, 100), p_Surv=p_surv)
# sim_flo <- data.frame(Flo=c(0:46), p_Flo=p_flo)
# sim_abund <- data.frame(Abund=c(0:70), p_Abund=p_abu)

ggplot(aes(x=Elev, y=p_Surv),data=sim_surv) + geom_smooth()
# ggplot(aes(x=Flo, y=p_Flo), data=sim_flo) + geom_point()
# ggplot(aes(x=Abund, y=p_Abund), data=sim_abund) + geom_point()





# #Posterior regression plots
# dev.off()
# source('post.plot.R')
# post.plot(Samples, param=c('beta.phi1.', 'beta.phi2.'))
# 
# # Bundle data
# dat <- read.table('./data/Solidago_PA_2008-2012.txt')
# dat %>% group_by(YEAR) %>% filter(occurrence>0) %>% summarize(dem=n())
# 
# grid <- expand.grid(long=dat$Easting,lat=dat$Northing)
# grid$occ <- rbinom(dim(grid)[1],1,prob = dat$occurrence) # occupancy
# grid$col <- ifelse(grid$occ==1,'orange','grey')
# 
# plot(lat~long,col=col,pch= 15, data=grid,ann=F,axes=F,cex=1.5)
# plot(lat~long,col='grey',pch= 15, data=grid,ann=F,axes=F,cex=1.5)
# 
# 
# Analyze_JAGSmodel <- function(x) {
#   # Plot the mcmc chain and the posterior sample for p
#   
#   Samples <- JAGS_model  
#   plot(Samples)
#   
#   # convergence check
#   gelman.diag(Samples)
#   
#   crosscorr.plot(Samples)
#   
#   # Statistical summaries of the posterior sample
#   summary(Samples)
#   
#   Samples.wd <- window(Samples, start = 1500 )
#   plot(Samples.wd)
#   
#   # Plot predicted quantities over a range of years
#   turnout.mat <- rbind(Samples.wd[[1]], Samples.wd[[2]], Samples.wd[[3]])
#   Samp.out <- data.frame(turnout.mat)
#   Samp.phis <- Samp.out[, grep("occ", colnames(Samp.out), fixed=T)]
#   
#   #load package
#   require(MCMCvis)
#   
#   MCMCplot(Samples, 
#            params = c('beta.phi1', 'beta.phi2', 'beta.phi3'),
#            horiz = FALSE,
#            ref_ovl = TRUE,
#            xlab = 'My x-axis label', 
#            main = 'MCMCvis plot',
#            labels=c('intercept', 'elevation', 'annual temp'))
#   
#   MCMCplot(Samples, 
#            params = c('a', 'y2'),
#            horiz = FALSE,
#            ref_ovl = FALSE,
#            xlab = 'My x-axis label', 
#            main = 'MCMCvis plot',
#            labels=c('distance decay', 'dispersal kernel'))
#   
#   #run summary function
#   Samp.sum <- MCMCsummary(Samples, round=3)
#   
#   ilogit <- function(u) { return(1 / (1 + exp(-u))); }
#   
#   #look at phis
#   library(ggplot2)
#   #OCC
#   occ<-Samp.sum[5:13690,1]
#   
#   
#   df_occ <- data.frame(cbind(ilogit(occ[1:2738]),
#                              ilogit(occ[2739:5476]), 
#                              ilogit(occ[5477:8214]),
#                              ilogit(occ[8215:10952]),
#                              ilogit(occ[10953:13690]),
#                              jags_data$elev[-rem.ind,1]))
#   colnames(df_occ)[1] <- "occ08";
#   colnames(df_occ)[2] <- "occ09";
#   colnames(df_occ)[3] <- "occ10";
#   colnames(df_occ)[4] <- "occ11";
#   colnames(df_occ)[5] <- "occ12";
#   colnames(df_occ)[6] <- "elev";
#   df <- melt(df_occ, id.vars = 'elev')
#   
#   dat <- data.frame(jags_data$C[-rem.ind,1:5])
#   dat<- cbind(dat, jags_data$elev[-rem.ind,1])
#   colnames(dat) <- c('Y1','Y2','Y3','Y4','Y5', 'elev')
#   dat <- melt(dat, id.vars = 'elev')
#   df$variable<-dat$variable
#   p1<-ggplot(dat, aes(x=elev, y=value)) +
#     geom_rug() + xlab('Elevation') + 
#     ylab('Occupancy probability')
#   p1 + geom_point(aes(df$elev, df$value), col='darkgrey')
#   # #C
#   # Col<-Samp.sum[1:5476,1]
#   # df_phi <- data.frame(cbind(Col[1:2738],Col[2739:5476], jags_data$cov1[-rem.ind,1]))
#   # colnames(df_phi)[1] <- "C1";
#   # colnames(df_phi)[2] <- "C2";
#   # colnames(df_phi)[3] <- "elev";
#   # df <- melt(df_phi, id.vars = 'elev')
#   
#   C_density_plot <-
#     ggplot(df, aes(x=value, col=variable)) +
#     geom_line(stat = "density", adjust=3) +
#     scale_x_continuous(name="probability of survival",
#                        limits=c(0,1), breaks=(0.2 * (0:5))) +
#     scale_y_continuous(name="probability density",
#                        limits=c(0,4.5), breaks=(0:4)) +
#     ggtitle("Posterior: Colonization (C)");
#   
#   plot(C_density_plot);
#   
#   ggplot(df, aes(x=elev, y=value, col=variable)) + geom_jitter()
#   geom_point() + geom_smooth()
#   
#   
#   
#   
#   #Plot survival based on blaaaaa
#   pred.DAY <- 1:102 # original dates span from Jan 18 - April 28, 2016; 102 days
#   pred.trDAY <- scale(pred.DAY)
#   modsims <- rstan::extract(OBFL.occ)
#   nsim <- length(modsims$lp__)  
#   newp <- array(dim=c(length(pred.DAY), nsim))
#   for(i in 1:nsim){ ## Loop through simulations 
#     surv[,i] <- plogis(Samples$beta.phi[i] + Samples$beta.phi1[i]*pred.trDAY + Samples$beta.phi2[i]*pred.trDAY^2)
#   }
#   
#   plot(NA, ylim=c(0,1), xlim=c(1, 102), axes=T, ylab="Detection probability", xlab="Julian Day (2016)", cex.lab=1.2)
#   polygon(c(pred.DAY, pred.DAY[length(pred.DAY):1]), c(apply(newp, 1, quantile, probs=0.025), apply(newp, 1, quantile, probs=0.975)[length(pred.DAY):1]), col="grey80", border="grey80")
#   lines(pred.DAY, apply(newp, 1, mean), col="blue", lwd=2.)
#   
#   
#   #EXTRAS  
#   ## Simulate the range of the moderating variable
#   
#   cov1.sim <- seq(min(cov1), max(cov1), by = 0.1)
#   
#   ## Calculate conditional effect of X1 across the range of X2
#   
#   ## Bayes:
#   
#   int.mcmc <- as.mcmc.list(Samples)
#   int.mcmc.mat <- as.matrix(int.mcmc)
#   int.sim <- matrix(rep(NA, nrow(int.mcmc.mat)*length(cov1.sim)), nrow = nrow(int.mcmc.mat))
#   for(i in 1:length(cov1.sim)){
#     int.sim[, i] <- int.mcmc.mat[2] + cov1.sim[i] + int.mcmc.mat[4] 
#   }
#   
#   ## Note: the variance now comes from the posterior, not the vcov matrix
#   
#   bayes.c.eff.mean <- apply(int.sim, 2, mean)
#   bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
#   bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))
#   
#   plot.dat <- data.frame(cov1.sim, bayes.c.eff.lower, bayes.c.eff.mean, bayes.c.eff.upper)
#   library(ggplot2)
#   
#   ## Use blue for Bayesian, red for frequentist estimates. Transparency to allow overlay; purple indicates complete overlay. Take a close look at the upper and lower limits of the CI for each estimate.
#   
#   ## Foundation for the plot & line for the posterior mean of the Bayesian conditional effect
#   p <- ggplot(plot.dat, aes(x = cov1.sim, y = bayes.c.eff.mean)) + geom_line(color = "blue", alpha = 0.8, size = 0.5)
#   
#   ## CI for the Bayesian conditional effect
#   p <- p + geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "blue", alpha = 0.2)
#   
#   ## Lines for the lower and upper bound of the Bayesian conditional effect
#   p <- p + geom_line(aes(x = cov1.sim, y = bayes.c.eff.lower), color = "blue", alpha = 0.8, size = 0.5) + 
#     geom_line(aes(x = cov1.sim, y = bayes.c.eff.upper), color = "blue", alpha = 0.8, size = 0.5)
#   
#   ## Plot labels and theme
#   p <- p + xlab("Altitude") + ylab("Conditional effect of X2") + theme_bw()
#   
#   ## Print the plot
#   p
#   
#   return(p) 
#   
# }
