#run.jags(module='glm')
# simplified spatial model

#### MODEL 1 ####
# This model removes most predictors and sets a fixed detection parameter at 0.95
sink("occ1.txt") 

cat(
  "

model{
  
  # Observation model
  for(t in 1:t.max){
    for(i in 1:n.sites){
       muY[i,t] <- z[i,t]*0.95
       y[i,t] ~ dbern(muY[i,t])
    }
  }
  
  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
    z[i,1] ~ dbern(psi1)
    
    for(t in 1:(t.max-1)){
      #  Model of local survival (1-extinction) at site i
      logit(phi[i,t]) <- beta_phi 
      C[i,t] <- beta_C
    
      Ez[i,t+1] = (1 - z[i,t])*C[i,t] + phi[i,t]*z[i,t] + 0.00001
      #For year 2, colonization if it was not occupied and survival if it was alive the previous year 
      #This won't normalize without the small additive constant. Why? 
      
      
      #new occupancy probability  
      z[i,t+1] ~ dbern(Ez[i,t+1])
    }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)
  #p ~ dunif(0, 1)
  
  #For predictors
  ##Survival
  beta_phi ~ dnorm(0,0.3)
  ##Colonization
  beta_C ~ dnorm(0,0.3)
  
  
  ###Derived quantities block
  psi[1]=psi1
  n_occ[1]=sum(z[1:n.sites,1])
  
  for(t in 1:(t.max-1)) {
    n_occ[t+1]=sum(z[1:n.sites,t+1])
    
	}
  }
  ",fill = TRUE, file="occ1.txt")

sink()

#Put in truncated dispersal kernel, estimate range limit


### 2) Set up a list that contains all the necessary data
Data_simple <- list(n.sites = nrow(occ), t.max = ncol(occ), y = occ)


# 3) Specify a function to generate inital values for the parameters
init.pa = occ
inits_fn = function() list(psi1 = 0.1, beta_phi=2, beta_C=2, z = init.pa)
#inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 2, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names = c("beta_phi","beta_C", "n_occ")  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)
    


#### MODEL 2 ####
# This model provides z as data with 1s fixed (now added randomly, need to get herb chrono data from gregoire). Then adds a p for detection probability, which I estimate. Is about 0.75!
sink("occ1.txt") 

cat(
  "
  
  model{
  
  # Observation model
  for(t in 1:t.max){
  for(i in 1:n.sites){
  muY[i,t] <- z[i,t]*p
  y[i,t] ~ dbern(muY[i,t])
  }
  }
  
  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)
  
  for(t in 1:(t.max-1)){
  #  Model of local survival (1-extinction) at site i
  logit(phi[i,t]) <- beta_phi 
  C[i,t] <- beta_C
  
  Ez[i,t+1] = (1 - z[i,t])*C[i,t] + phi[i,t]*z[i,t] + 0.00001
  #For year 2, colonization if it was not occupied and survival if it was alive the previous year 
  #This won't normalize without the small additive constant. Why? 
  
  
  #new occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)
  p ~ dunif(0, 1)
  
  #For predictors
  ##Survival
  beta_phi ~ dnorm(0,0.3)
  ##Colonization
  beta_C ~ dnorm(0,0.3)
  
  
  ###Derived quantities block
  psi[1]=psi1
  n_occ[1]=sum(z[1:n.sites,1])
  
  for(t in 1:(t.max-1)) {
  n_occ[t+1]=sum(z[1:n.sites,t+1])
  
  }
  }
  ",fill = TRUE, file="occ1.txt")

sink()

#Put in truncated dispersal kernel, estimate range limit


### 2) Set up a list that contains all the necessary data
z<-occ
z[z==0]<-NA
z[sample(which(rowSums(z, na.rm=T)>=1), 50),]<-1

Data_simple <- list(n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z)


# 3) Specify a function to generate inital values for the parameters
#init.pa = occ
inits_fn = function() list(psi1 = 0.1, beta_phi=2, beta_C=2, z = z, p = 0.9)
#inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 2, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names = c("beta_phi","beta_C", "n_occ", 'p')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)




#### MODEL 3 ####

#Include p and a false postiive test with priors (stronger prior on the pn)
sink("occ1.txt") 

cat(
  "
  
  model{
  
  # Observation model
  for(t in 1:t.max){
  for(i in 1:n.sites){
  muY[i,t] <- (z[i,t]*p) + (1-z[i,t])*pn
  y[i,t] ~ dbern(muY[i,t])
  }
  }
  
  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)
  
  for(t in 1:(t.max-1)){
  #  Model of local survival (1-extinction) at site i
  logit(phi[i,t]) <- beta_phi 
  C[i,t] <- beta_C
  
  Ez[i,t+1] = (1 - z[i,t])*C[i,t] + phi[i,t]*z[i,t] + 0.00001
  #For year 2, colonization if it was not occupied and survival if it was alive the previous year 
  #This won't normalize without the small additive constant. Why? 
  
  
  #new occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)
  p ~ dbeta(1, 1)
  pn ~ dbeta(20, 1)
  
  #For predictors
  ##Survival
  beta_phi ~ dnorm(0,0.3)
  ##Colonization
  beta_C ~ dnorm(0,0.3)
  
  
  ###Derived quantities block
  psi[1]=psi1
  n_occ[1]=sum(z[1:n.sites,1])
  
  for(t in 1:(t.max-1)) {
  n_occ[t+1]=sum(z[1:n.sites,t+1])
  
  }
  }
  ",fill = TRUE, file="occ1.txt")

sink()

#Put in truncated dispersal kernel, estimate range limit


### 2) Set up a list that contains all the necessary data
z<-occ
z[z==0]<-NA
z[sample(which(rowSums(z, na.rm=T)>=1), 50),]<-1

Data_simple <- list(n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z)


# 3) Specify a function to generate inital values for the parameters
#init.pa = occ
inits_fn = function() list(psi1 = 0.1, beta_phi=2, beta_C=2, z = z, p = 0.5, pn=0.9)
#inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 2, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names = c("beta_phi","beta_C", "n_occ", 'p', 'pn')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)



#### MODEL 4 ####

#Include missing data in predictors (test!)
sink("occ1.txt") 

cat(
  "
  
  model{
  
  # Observation model
  for(t in 1:t.max){
  for(i in 1:n.sites){
  muY[i,t] <- (z[i,t]*p) + (1-z[i,t])*pn
  y[i,t] ~ dbern(muY[i,t])
  }
  }
  
  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)


  #Model of missing data for elev
  elev[i] ~ dnorm(muelev, telev) #outside t loop so defined only once on left hand side
  
  for(t in 1:(t.max-1)){

  #  Model of local survival (1-extinction) at site i
  logit(phi[i,t]) <- beta_phi1 + beta_phi2*elev[i] 
  C[i,t] <- beta_C
  
  Ez[i,t+1] = (1 - z[i,t])*C[i,t] + phi[i,t]*z[i,t] + 0.00001
  #For year 2, colonization if it was not occupied and survival if it was alive the previous year 
  #This won't normalize without the small additive constant. Why? 
  
  
  #new occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)
  p ~ dbeta(1, 1)
  pn ~ dbeta(20, 1)
  
  #For missing elev data
  muelev ~ dnorm(0, 10^-6)
  telev ~ dgamma(10,10)
  
  #For predictors
  ##Survival
  beta_phi1 ~ dnorm(0,0.3)
  beta_phi2 ~ dnorm(0,0.3)
  ##Colonization
  beta_C ~ dnorm(0,0.3)
  
  
  ###Derived quantities block
  psi[1]=psi1
  n_occ[1]=sum(z[1:n.sites,1])
  
  for(t in 1:(t.max-1)) {
  n_occ[t+1]=sum(z[1:n.sites,t+1])
  
  }
  }
  ",fill = TRUE, file="occ1.txt")

sink()

#Put in truncated dispersal kernel, estimate range limit


### 2) Set up a list that contains all the necessary data
z<-occ
z[z==0]<-NA
z[sample(which(rowSums(z, na.rm=T)>=1), 50),]<-1

#add some missing data
newelev<- elev
newelev[sample(1:2000, 50),sample(1:4, 2)] <- NA 
newelev<- newelev[,1]

Data_simple <- list(n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z, elev=newelev)


# 3) Specify a function to generate inital values for the parameters
#init.pa = occ
inits_fn = function() list(psi1 = 0.1, muelev = 1000, telev= 10^-3, beta_phi1=2, beta_phi2=2, beta_C=2, z = z, p = 0.5, pn=0.9)
#inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 2, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names = c("beta_phi1","beta_phi2", "beta_C", "n_occ", 'p', 'pn')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)


#### MODEL 5 ####

sink("occ1.txt") 
#Includes correct model (full model) with missing data in dem and flo, but no dispersal kernel

cat(
  "
  
  model{
  
  # Observation model
  for(t in 1:t.max){
  for(i in 1:n.sites){
  muY[i,t] <- z[i,t]*p
  y[i,t] ~ dbern(muY[i,t])
  }
  }
  
  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)

    
    for(t in 1:(t.max-1)){
  
    #Model of missing data for dem and flo (in t loop because needs to be defined multiple times
    flo[i,t] ~ dnorm(mu_flo, tau_flo) 
    dem[i,t] ~ dnorm(mu_dem, tau_dem) 

    #  Model of local survival (1-extinction) at site i
    logit(phi[i,t]) <- beta_phi1 + beta_phi2*elev[i] 
    
    # Colonization prob
    gamma[i,t] = beta_gamma #pairs of sites, integrate over all sites.
  
    # Flowering:
    logit(kappa[i,t]) = beta_f[1] + beta_f[2] * flo[i,t] 
  
    # Relative abundance
    logit(lambda[i,t]) = beta_a[1] + beta_a[2] * dem[i,t]
  
    Ez[i,t+1] = kappa[i,t]*lambda[i,t]*gamma[i,t]*(1-z[i,t]) + phi[i,t]*z[i,t]+0.0001
    #For year 2, colonization if it was not occupied and survival if it was alive the previous year 
  
    #new occupancy probability  
    z[i,t+1] ~ dbern(Ez[i,t+1])
    }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)
  p ~ dbeta(1, 1)
  
  #For missing data in flowering and dem data
  mu_flo ~ dgamma(10, 10) #won't be negative! (maybe use a total vs. average <- poisson for the respoonse)
  tau_flo ~ dgamma(10,10) #so could use the total counts, use dgamma as conjugate to poisson for mean
  mu_dem ~ dgamma(10, 10)
  tau_dem ~ dgamma(10,10)
  

  #For predictors
  ##Survival
  beta_phi1 ~ dnorm(0,0.3)
  beta_phi2 ~ dnorm(0,0.3)
  ##Colonization
  beta_gamma ~ dnorm(0,0.3)
  
  ###Flowering
  beta_f[1] ~ dnorm(0,0.3)
  beta_f[2] ~ dnorm(0, 0.3)
  
  ###Relative abundance
  beta_a[1] ~ dnorm(0,0.3)
  beta_a[2] ~ dnorm(0, 0.3)
  
  ###Derived quantities block
  psi[1]=psi1
  n_occ[1]=sum(z[1:n.sites,1])
  
  for(t in 1:(t.max-1)) {
  n_occ[t+1]=sum(z[1:n.sites,t+1])
  
  }
  }
  ",fill = TRUE, file="occ1.txt")

sink()


### 2) Set up a list that contains all the necessary data
z<-occ
z[z==0]<-NA
z[sample(which(rowSums(z, na.rm=T)>=1), 50),]<-1

Data_simple <- list(n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z, dem=dem, flo=flo, elev=elev[,1])

# 3) Specify a function to generate inital values for the parameters
#init.pa = occ
inits_fn = function() list(psi1 = 0.1, mu_dem = 10, tau_dem= 10^-3, mu_flo = 10, tau_flo= 10^-3, beta_phi1=2, beta_phi2=2, beta_gamma=2, z = z, p = 0.5)
#inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 2, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names = c("n_occ", 'p')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)

p_values <- c(Samples[[1]][,'p'], Samples[[2]][,'p']) #[[], grep("p", colnames(Samples), fixed=T)]]
median(p_values)




#### MODEL 6 ####

sink("occ1.txt") 
#Includes correct model (full model) with missing data in dem and flo, but no dispersal kernel, plus site and time random effect!

cat(
  "
  
  model{
  
  # Observation model
  for(t in 1:t.max){
  for(i in 1:n.sites){
  muY[i,t] <- z[i,t]*p
  y[i,t] ~ dbern(muY[i,t])
  }
  }
  
  for(t in 1:(t.max-1)){
  trand[t] ~ dnorm... 0, taut} #if a time point is bad, it's bad for all sites so keep additive
  
  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)
  
  rand ~ dnorm(0, tau)
  for(t in 1:(t.max-1)){
  
  #Model of missing data for dem and flo (in t loop because needs to be defined multiple times
  flo[i,t] ~ dnorm(mu_flo, tau_flo) #dpois if total, may likely need a negative binomial (search help for dnbin expects number of failures before successes, need to convert it). Or if you have enough observations, you could assume this is normal and the variances are estimated independently (vs. the poisson)
  dem[i,t] ~ dnorm(mu_dem, tau_dem) 
  
  #  Model of local survival (1-extinction) at site i
  logit(phi[i,t]) <- beta_phi1 + beta_phi2*elev[i] + rand[i] + trand[t] #if you wanted to simulate for a new time point, would be adding some noise to survival probability. You can say time effect is 0 for new time points to assume it's average.But makes more sense to just let JAGs do it, gives more uncertainty in following year but that makes sense. Makes sense for sites too.
  
  # Colonization prob
  gamma[i,t] = beta_gamma #pairs of sites, integrate over all sites.
  
  # Flowering:
  logit(kappa[i,t]) = beta_f[1] + beta_f[2] * flo[i,t] 
  
  # Relative abundance
  logit(lambda[i,t]) = beta_a[1] + beta_a[2] * dem[i,t]
  
  Ez[i,t+1] = kappa[i,t]*lambda[i,t]*gamma[i,t]*(1-z[i,t]) + phi[i,t]*z[i,t]+0.0001
  #For year 2, colonization if it was not occupied and survival if it was alive the previous year 
  
  #new occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)
  p ~ dbeta(1, 1)
  
  #For missing data in flowering and dem data
  mu_flo ~ dgamma(10, 10) #won't be negative! (maybe use a total vs. average <- poisson for the respoonse)
  tau_flo ~ dgamma(10,10) #so could use the total counts, use dgamma as conjugate to poisson for mean
  mu_dem ~ dgamma(10, 10)
  tau_dem ~ dgamma(10,10)
  
  #conugate for mean paramater for normal dist is dnorm, for the precision its dgamma.
  #For predictors
  ##Survival
  beta_phi1 ~ dnorm(0,0.3)
  beta_phi2 ~ dnorm(0,0.3)
  ##Colonization
  beta_gamma ~ dnorm(0,0.3)
  
  ###Flowering
  beta_f[1] ~ dnorm(0,0.3)
  beta_f[2] ~ dnorm(0, 0.3)
  
  ###Relative abundance
  beta_a[1] ~ dnorm(0,0.3)
  beta_a[2] ~ dnorm(0, 0.3)
  
  ###Derived quantities block
  psi[1]=psi1
  n_occ[1]=sum(z[1:n.sites,1])
  
  for(t in 1:(t.max-1)) {
  n_occ[t+1]=sum(z[1:n.sites,t+1])
  
  }
  }
  ",fill = TRUE, file="occ1.txt")

sink()
#you could ask just for the mean, is in rjags! Look at means only, not posterior (example in WAIC folder in Day 3)
#Put in truncated dispersal kernel, estimate range limit


### 2) Set up a list that contains all the necessary data
z<-occ
z[z==0]<-NA
z[sample(which(rowSums(z, na.rm=T)>=1), 50),]<-1

Data_simple <- list(n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z, dem=dem, flo=flo, elev=elev[,1])

# 3) Specify a function to generate inital values for the parameters
#init.pa = occ
inits_fn = function() list(psi1 = 0.1, mu_dem = 10, tau_dem= 10^-3, mu_flo = 10, tau_flo= 10^-3, beta_phi1=2, beta_phi2=2, beta_gamma=2, z = z, p = 0.5)
#inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 2, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names = c("n_occ", 'p')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)

