#Individual Patch Performance Models
#Using dem_data from readdata_jags_dem
library('runjags')

head(dem_data)


#### Density over time ####
N<-nrow(dem_data)
center_apply <- function(y) y - mean(y)
X<-center_apply(elev[,1])
Y<-as.matrix(dem_data[,c('Density08', 'Density10', 'Density12')])
dat <- list(N=N, X=X, Y=Y[,1])
  
cat('model {
    intercept ~ dnorm(0, 1/10^6) #remember sd is precision
    coef ~ dnorm(0, 1/10^6)
    tau ~ dgamma(10,10)

    for(i in 1:N){
      Y[i] ~ dnorm(mu[i], tau)
      mu[i] <- intercept + coef*X[i] 
    }

    #data# Y, X, N
    #monitor# intercept, coef, tau
    }', file='modelfile.txt')

# a named list or a function

results <- run.jags(model='modelfile.txt', data = dat, burnin=2000, sample=10000, n.chains=2)

plot(results, 'trace')

summary(results)			


Y<-as.matrix(dem_data[,c('Density08', 'Density10', 'Density12')])
dat <- list(N=N, X=X, Y=Y)


cat('model {
    intercept ~ dnorm(0, 1/10^6) #remember sd is precision
    coef ~ dnorm(0, 1/10^6)
    tau ~ dgamma(0.001,0.001)
    ttau ~ dgamma(0.001,0.001)

    for(t in 1:3){
      trand[t] ~ dnorm(0, ttau)
    }
    
    for(i in 1:N){
        for(t in 1:3) {
          Y[i,t] ~ dnorm(mu[i,t], tau)
          log(mu[i,t]) <- intercept + coef*X[i]  + trand[t]
          Y_rep[i,t] ~ dnorm(mu[i,t], tau)
        }
    }

    #data# Y, X, N
    #monitor# intercept, coef, tau, ttau, Y_rep
    }', file='modelfile.txt')


results <- run.jags(model='modelfile.txt', data = dat, burnin=2000, sample=2000, thin=2, n.chains=2)

plot(results)

library(bayesplot)
library(Bayesfactor)
posterior(results) %>% 
  gf_point(intercept ~ coef, color = ~ chain, alpha = 0.2, size = 0.4) %>%
  gf_density2d(alpha = 0.5)

#Extract samples
coef <- combine.mcmc(results, var='coef')
intercept <- combine.mcmc(results, var='intercept')

#Prepare samples for plotting
df_obs <- data.frame(k1=Y[,1], k2=Y[,2], k3=Y[,3], elev=X)

df_pred <- data.frame(k1 = Y_rep[,1],
                      k2 = Y_rep[,2],
                      k3 = Y_rep[,3],
                      elev=X)

m <- melt(table(df_pred), value.name = "dem")

#Then plot
p1 <- ggplot(m, aes(x = elev, y = k3)) +
  geom_point(shape=19,
             size = m[,'dem']/5,
             alpha=0.8,
             colour = "grey") +
  geom_point(data = df_obs, aes(x = elev,
                                y = k3),
             size = 6,
             shape = 19,
             colour = "red") +
  xlab("Success count 1") +
  ylab("Success count 2") +
  theme_bw(base_size = 12, base_family = "Helvetica")

#Or other tutorial plot
plotModelOutput = function(jagsmodel, Y) {
  # attach the model
  attach.jags(results)
  x = seq(1,length(Y))
  summaryPredictions = cbind(apply(predY,2,quantile,0.025), apply(predY,2,mean),
                             apply(predY,2,quantile,0.975))
  plot(Y, col="white",ylim=c(min(c(Y,summaryPredictions)),max(c(Y,summaryPredictions))),
       xlab="",ylab="95% CIs of predictions and data",main=paste("JAGS results:",
                                                                 jagsmodel$model.file))
  polygon(c(x,rev(x)), c(summaryPredictions[,1], rev(summaryPredictions[,3])),
          col="grey70",border=NA)
  lines(summaryPredictions[,2])
  points(Y)
}
and we can use the function to plot the predicted posterior mean with 95%
CIs, as well as the raw data. For example, try
plotModelOutput(mod_lmcor_intercept, Wind)

#### FLOWERS (over time?) ####
Y<-as.matrix(dem_data[,c('Flower10', 'Flower12')])
dat <- list(N=N, X=X, Y=Y)

cat('model {
    intercept ~ dnorm(0, 1/10^6) #remember sd is precision
    coef ~ dnorm(0, 1/10^6)
    tau ~ dgamma(10,10)
    
    #ntau~ dgamma(10,10)
    #ttau ~ dgamma(10,10)
    
    #for(t in 1:(3)){
    #trand[t] ~ dnorm(0, ttau)
    #}
    
    for(i in 1:N){
    #rand ~ dnorm(0, ntau)
    for(t in 1:2) {
    Y[i,t] ~ dnorm(mu[i,t], tau)
    mu[i,t] <- intercept + coef*X[i] #+ rand[i] + trand[t]
    }
    }
    #data# Y, X, N
    #monitor# intercept, coef, tau
    }', file='modelfile.txt')


results <- run.jags(model='modelfile.txt', data = dat, burnin=2000, sample=10000, n.chains=2)

plot(results)

results			

coef <- combine.mcmc(results, var='coef')
intercept <- combine.mcmc(results, var='intercept')

y<- intercept + coef*X

plot(Y[,2] ~ X)


#### Density (over time) JAGS template ####
center_apply <- function(y) y - mean(y)
dat <- data.frame(dem_data[,c('Density08', 'Density10', 'Density12')])
elev<-center_apply(Data_simple$elev)
dat$elev <- elev
dat$patch <- c(1:2741)
dat <- dat %>% gather(year, dem, -elev, -patch)
model <- template.jags(dem ~ elev + (1 | year), data=dat, family='gaussian')
results <- run.jags("JAGSmodel.txt", method="int", sample = 5000, modules = 'glm')
results

plot(results)

# Having #residual# and #fitted# in the JAGS specification allows for:
model.residuals <- residuals(results, output='mean')
model.fitted <- fitted(results, output='mean')
plot(model.fitted, data$dem)
plot(model.fitted, model.residuals); abline(h=0)

#### EXTRA CODE ####

#Including time as a random effect, but number of observations < 5 so NOT A GOOD IDEA DUE TO MASSIVE SHRINKAGE!!! 

# dat <- list(N=N, X=X, Y=Y)
# cat('model {
#     intercept ~ dnorm(0, 1/10^6) #remember sd is precision
#     coef ~ dnorm(0, 1/10^6)
#     tau ~ dgamma(10,10)
#     
#     #ntau~ dgamma(10,10)
#     #ttau ~ dgamma(10,10)
#     
#     #for(t in 1:(3)){
#     #trand[t] ~ dnorm(0, ttau)
#     #}
#     
#     for(i in 1:N){
#     #rand ~ dnorm(0, ntau)
#     for(t in 1:3) {
#     Y[i,t] ~ dnorm(mu[i,t], tau)
#     log(mu[i,t]) <- intercept + coef*X[i] #+ rand[i] + trand[t]
#     }
#     }
#     #data# Y, X, N
#     #monitor# intercept, coef, tau
#     }', file='modelfile.txt')
# 
# 
# results <- run.jags(model='modelfile.txt', data = dat, burnin=2000, sample=10000, n.chains=2)
# 
