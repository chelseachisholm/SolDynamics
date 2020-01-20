#### Process Output from JAGS models ####
#07.06.2019 C. Chisholm

#load(occ,od_21122018.RData)

source('./code/helperfunctions.R')

### Convergence checks ####
# Plot the mcmc chain and the posterior sample for p
plot(Samples, ask=T)

# convergence check
gelman.plot(Samples)
coda::gelman.diag(Samples)

# Statistical summaries of the posterior sample
summary(Samples)

Samples.wd <- window(Samples, start = 1500 )
plot(Samples.wd)

#### Plot 1: Change in number of patches ####
#mcmc <- as.data.frame(cbind(Samples[[1]], Samples[[2]], Samples[[3]])) 
mcmc <- as.data.frame(cbind(Samples[[1]])) 
coefs = mcmc[, c("n_occ[1]", "n_occ[2]", "n_occ[3]", "n_occ[4]", "n_occ[5]")]

#For now the survival and number of occupied patches is way too high. What's going on?

#### Plot 1: Survival and Colonisation ####

nelev<- scale(elev)
mcmc <- as.data.frame(cbind(Samples[[1]], Samples[[2]], Samples[[3]])) 

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

ggplot(dat, aes(y = estimate, x = elev)) + geom_line() + 
  geom_ribbon(aes(ymin = l_cred, ymax = u_cred), fill = "blue", alpha = 0.3) + 
  scale_y_continuous("Survival") +
  scale_x_continuous("Elevation (m)") + theme_classic()

ggplot(dat, aes(y = estimate, x = elev)) + geom_line() + 
  scale_y_continuous("Survival") +
  scale_x_continuous("Elevation") + theme_classic(base_size = 16)



#Using ggplot and matrix mulitplication
newdata = data.frame(x = seq(min(nelev, na.rm = TRUE), max(nelev, na.rm = TRUE),
                             len = 100))
Xmat = model.matrix(~x, newdata)
coefs = mcmc[, c("beta_phi1", "beta_phi2")]
fit = antilogit(coefs * Xmat)
newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "HPDinterval"))
ggplot(newdata, aes(y = estimate, x = x)) + geom_line() + geom_ribbon(aes(ymin = conf.low,
                                                                          ymax = conf.high), fill = "blue", alpha = 0.3) + scale_y_continuous("Y") +
  scale_x_continuous("X") + theme_classic()


#### Plot 2: Demography Effects ####
#logit(phi[i,t]) <- beta_phi[1] + beta_phi[2] * elev[i,t]
# Flowering:
#logit(kappa[i,t]) = beta_f[1] + beta_f[2] * flo[i,t] 
# Relative abundance
#logit(lambda[i,t]) = beta_a[1] + beta_a[2] * dem[i,t]

p_surv = antilogit(par_vals[5,1] + par_vals[6,1]*elev)
p_flo = antilogit(par_vals[3,1] + par_vals[4,1]*flo)
p_abu = antilogit(par_vals[1,1] + par_vals[2,1]*dem)
Covdat <-data.frame(Year=as.factor(rep(1:5, each=2741)), flo=c(flo), elev=c(elev), dem=c(dem), Surv=c(p_surv), Flo=c(p_flo), Abund=c(p_abu))

ggplot(aes(x=elev, y=Surv), data=Covdat) + geom_smooth() + geom_line() + facet_wrap(~Year)
ggplot(aes(x=elev, y=flo), data=Covdat) + geom_point() #+ geom_line() #+ facet_wrap(~Year)
ggplot(aes(x=dem, y=Abund), data=Covdat) + geom_point() + geom_line() + facet_wrap(~Year)

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

#### Plot landscape connectivity ####

par_vals = summary(Samples)$statistics
gammas <- par_vals %>% as.data.frame() %>% rownames_to_column(var='param') %>%
  filter(grepl('gamma', param)) %>% filter(param != 'gamma0') %>%
  mutate(Year= rep(1:4, each=2741), Patch=rep(1:2741, times=4), Elev=rep(elev[,1], times=4))

gammas %>% group_by(Year) %>% summarize(mean=mean(Mean))
ggplot(gammas, aes(x=Elev, y=Mean)) + geom_point() + facet_wrap(~Year) 

ggplot(dat, aes(y = estimate, x = elev)) + geom_line() + 
  scale_y_continuous("Survival") +
  scale_x_continuous("Elevation") + theme_classic(base_size = 16) +
  geom_point(data=gammas, aes(x=Elev, y=Mean), col='grey')

gamma0 <- par_vals %>% as.data.frame() %>% rownames_to_column(var='param') %>%
  filter(grepl('gamma', param)) %>% filter(param == 'gamma0') 


#### Plot predicted quantities with uncertainty ####

surv=NULL
nsim=1000
  for(i in 1:nsim){ ## Loop through simulations
    surv[i] <- antilogit(X$beta_phi1[i] + X$beta_phi2[i]*elev[1])
  }

plot(NA, ylim=c(0,1), xlim=c(483:2052), axes=T, ylab="Survival probability", xlab="Elevation (m)", cex.lab=1.2)
polygon(c(surv, surv[length(surv):1]), c(apply(surv, 1, quantile, probs=0.025), apply(surv, 1, quantile, probs=0.975)[length(surv):1]), col="grey80", border="grey80")
lines(elev, apply(surv, 1, mean), col="blue", lwd=2.)
#par(mfrow = c(2,2))

#### Plot paramters per year? ####
X <- as.data.frame(cbind(Samples[[1]], Samples[[2]])) #beta_phi1 and 2
means <- apply(X,2,mean)
upperCI <- apply(X,2,quantile,0.975)
lowerCI <- apply(X,2,quantile,0.025)


nYears <- ncol(occ)
for(i in 1:length(means)) {
  plot(means[i],lwd=3,ylim=range(c(lowerCI[i],upperCI[i])),
       type="n")
  polygon(c(1:nYears),
          c(upperCI[i],rev(lowerCI[i]),upperCI[i]),col="skyblue",lty=0)
  lines(means[i],lwd=3)
}

detach.jags()



#Plot predicted quantities over a range of years?
turnout.mat <- rbind(Samples[[1]], Samples[[2]], Samples[[3]])
Samp.out <- data.frame(turnout.mat)
p <- Samp.out$beta.phi2

p.mean <- mean(p)

plot(p.mean ~ c(1:5), xlim = c(1,5))

list(int90 = HPDinterval(Samples[[1]], prob=.90),
     int95= HPDinterval(Samples[[1]], prob=.95),
     int99= HPDinterval(Samples[[1]], prob=.99))

dic.samples(jagsModel, 1000)

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