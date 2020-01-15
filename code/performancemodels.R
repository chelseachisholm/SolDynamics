#Individual Patch Performance Models
#Using dem_data from readdata_jags_dem
library('runjags')

head(dem_data)


#### Density over time ####
N<-nrow(dem_data)
center_apply <- function(y) y - mean(y)
X<-center_apply(dem_data$Altitude)
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
dat <- list(N=N, X=X, Y=Y[,1])


plot(results)

results			

cat('model {
    intercept ~ dnorm(0, 1/10^6) #remember sd is precision
    coef ~ dnorm(0, 1/10^6)
    tau ~ dgamma(10,10)

    #ntau~ dgamma(10,10)
    #ttau ~ dgamma(10,10)

    #for(t in 1:3){
    #trand[t] ~ dnorm(0, ttau)
    #}
    
    for(i in 1:N){
    #rand[i] ~ dnorm(0, ntau)
        for(t in 1:3) {
          Y[i,t] ~ dnorm(mu[i,t], tau)
          log(mu[i,t]) <- intercept + coef*X[i] #+ rand[i] #+ trand[t]
        }
    }
    #data# Y, X, N
    #monitor# intercept, coef, tau
    }', file='modelfile.txt')


results <- run.jags(model='modelfile.txt', data = dat, burnin=2000, sample=20000, n.chains=2)

plot(results)

coef <- combine.mcmc(results, var='coef')
intercept <- combine.mcmc(results, var='intercept')

y<- exp(intercept + coef*X)

plot(Y[,3] ~ X)

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
coda::HPDinterval(coef)

# Based on the full posteriors:
expectation <- intercept + coef * mean(X)
median(expectation)
coda::HPDinterval(expectation)


library(tidybayes)
library(modelr)
plot(Y[,1] ~ X)
lambda<-median(intercept)+median(coef)*X
lines(X, lambda)
library("ggplot2")
eq<- function(x) {median(intercept) + median(coef)*x}
dfpred <- data.frame(x=c(-300,500))
ggplot(dfpred, aes(x=x)) + stat_function(fun=eq, geom="line") + xlab("x") + ylab("y") +
  geom_point(data = as.data.frame(dat), aes(x=X, y=Y[,2]),size = 2) + 
  scale_x_continuous(c(-300,500)) 

rar <- as.data.frame(combine.mcmc(results)) %>% gather(key='parameter')
rar %>%
  ggplot(aes(x = value, y=as.factor(parameter))) +
  tidybayes::stat_intervalh()

rar %>%
  ggplot(dfpred, aes(x = x)) +
  stat_function(fun=eq, geom="line") + xlab("x") + ylab("y") +
  stat_lineribbon(fun=eq, .width = .95, color = "#08519C") +
  geom_point(data = rar, size = 2) +
  scale_fill_brewer()

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
