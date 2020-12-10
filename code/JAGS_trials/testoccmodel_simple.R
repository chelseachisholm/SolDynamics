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
# Detection probability!
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
#Use z from readdata_jags_dem.R, below is just a sample I tested
#For simulated data
# z<-occ
# z[z==0]<-NA
# z[sample(which(rowSums(z, na.rm=T)>=1), 50),]<-1


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
  muY[i,t] <- (z[i,t]*p) + (1-z[i,t])*pn #false negative p and false positive pn (stupid name)
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
#Use z from model 2 or readdata_jags_dem.R

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
  
  #Random year effect for each year
  for(t in 1:(t.max-1)){
  trand[t] ~ dnorm(0, taut) } #if a time point is bad, it's bad for all sites so keep additive
  
  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)
  
  #Random site effect
  rand[i] ~ dnorm(0, tau)

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

  #For random effects of year and time
  taut ~ dgamma(0.1,0.1)
  tau ~ dgamma(0.1,0.1)
  
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

Data_simple <- list(n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z, dem=dem, flo=flo, elev=elev[,1])

# 3) Specify a function to generate inital values for the parameters
#init.pa = occ
inits_fn = function() list(psi1 = 0.1, mu_dem = 10, tau_dem= 10^-3, mu_flo = 10, tau_flo= 10^-3, beta_phi1=2, beta_phi2=2, beta_gamma=2, z = z, p = 0.5)
#inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 2, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names = c("n_occ", 'p', 'tau', 'taut')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)

#### MODEL 7 ####

#Includes simplified model 1 with dispersal kernel!
#Two ways to do distance matrix. Simply exponential decay or a defined 'source strength'

#Create list of neighbours within max distance (using 500 for now)
#Create a list of neighbours within max.dist for each patch
max.dist <- 50
NB.list <- apply(distmat,1, function(x) which(x <= max.dist))
# Number of neighbours per patch
n.nb <- sapply(NB.list,length)
summary(n.nb)
# 
# Here I create a matrix NB.mat with the patch numbers (indices) of all neighbours
# for each patch
# and a re-formated distance matrix D.nb with the distances to these neighbours
NB.mat <- matrix(NA, nrow = nrow(occ), ncol = max(n.nb))
D.nb <- NB.mat
for(i in 1:nrow(occ)){
  NB.mat[i,1:n.nb[i]] <- NB.list[[i]]
  D.nb[i,1:n.nb[i]] <- distmat[i,NB.list[[i]]]
}

#  Pairwise 'source strength' calc and colonisation probability
# for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
#  gammaDistPairs[i,n,t] <- z[NB.mat[i,n],t] * alpha *exp(-alpha * D.nb[i, n])
# }
# 
# gamma[i,t] = 1 - prod(1-gammaDistPairs[i,,t])
# S[i,t] <- sum(S.nb[i,t,1:n.nb[i]])

#Or full pairwise model
#Pairwise colonization probability based on dispersal kernel 
# for(n in 1:nNeighbours){ 
#   #Exponential kernel, probability of dispersal to a site based on distance
#   gammaDistPairs[i,n,t] = gamma0 * exp(-gamma0*dist.all[i,n]) * z[n,t] 
# }#add in repro ability for other patches
# 
# # Colonization prob
# gamma[i,t] = 1 - prod(1-gammaDistPairs[i,,t])

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
  
  #  Pairwise 'source strength' calc and colonisation probability
  for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
  gammaDistPairs[i,n,t] <- 1-gamma0*exp(-gamma0*D.nb[i,n])*z[NB.mat[i,n],t] #This indexes the neighbour in question (NB.mat[i,n])
  }
  
  # Colonization prob
  gamma[i,t] = 1 - prod(1-gammaDistPairs[i,1:n.nb[i],t])
  
  
  ## Pr(z=1|z=1) = 1 - Pr(extinction)*Pr(not colonized) = 1 - ((1-phi)*(1-gamma))
  Ez[i,t+1] = (1 - z[i,t])*gamma[i,t] + (1-(1-phi[i,t])*(1-gamma[i,t]))*z[i,t] 
  #For year 2, colonization if it was not occupied and survival if it was alive the previous year 
  
  #new occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)

  gamma0 ~ dbeta(1,1)

  
  #For predictors
  ##Survival
  beta_phi ~ dnorm(0,0.3)
  ##Colonization
  #beta_C ~ dnorm(0,0.3)
  
  
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
Data_simple <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(occ), t.max = ncol(occ), y = occ)


# 3) Specify a function to generate inital values for the parameters
init.pa = occ
inits_fn = function() list(gamma0=2, psi1 = 0.1, beta_phi=2, z = z)
inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 2, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names = c("p", "beta_phi","gamma0", "n_occ")  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples)


#### MODEL 8 ####

#Includes correct model (full model) with missing data in dem and flo, plus site and time random effect and dispersal kernel!
#Two ways to do distance matrix. Simply exponential decay or a defined 'source strength'

#Create list of neighbours within max distance (using 500 for now)
#Create a list of neighbours within max.dist for each patch
max.dist <- 200
NB.list <- apply(distmat,1, function(x) which(x <= max.dist))
# Number of neighbours per patch
n.nb <- sapply(NB.list,length)
summary(n.nb)
# 
# Here I create a matrix NB.mat with the patch numbers (indices) of all neighbours
# for each patch
# and a re-formated distance matrix D.nb with the distances to these neighbours
NB.mat <- matrix(NA, nrow = nrow(occ), ncol = max(n.nb))
D.nb <- NB.mat
for(i in 1:nrow(occ)){
  NB.mat[i,1:n.nb[i]] <- NB.list[[i]]
  D.nb[i,1:n.nb[i]] <- distmat[i,NB.list[[i]]]
}

#  Pairwise 'source strength' calc and colonisation probability
# for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
#  gammaDistPairs[i,n,t] <- z[NB.mat[i,n],t] * alpha *exp(-alpha * D.nb[i, n])
# }
# 
# gamma[i,t] = 1 - prod(1-gammaDistPairs[i,,t])
# S[i,t] <- sum(S.nb[i,t,1:n.nb[i]])

#Or full pairwise model
#Pairwise colonization probability based on dispersal kernel 
# for(n in 1:nNeighbours){ 
#   #Exponential kernel, probability of dispersal to a site based on distance
#   gammaDistPairs[i,n,t] = gamma0 * exp(-gamma0*dist.all[i,n]) * z[n,t] 
# }#add in repro ability for other patches
# 
# # Colonization prob
# gamma[i,t] = 1 - prod(1-gammaDistPairs[i,,t])

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
  
  for(t in 1:(t.max-1)){
  trand[t] ~ dnorm(0, taut) } #if a time point is bad, it's bad for all sites so keep additive
  
  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)
  
  rand[i] ~ dnorm(0, tau)
  
  for(t in 1:(t.max-1)){
  
  #Model of missing data for dem and flo (in t loop because needs to be defined multiple times
  flo[i,t] ~ dnorm(mu_flo, tau_flo) #dpois if total, may likely need a negative binomial (search help for dnbin expects number of failures before successes, need to convert it). Or if you have enough observations, you could assume this is normal and the variances are estimated independently (vs. the poisson)
  dem[i,t] ~ dnorm(mu_dem, tau_dem) 
  
  #  Model of local survival (1-extinction) at site i
  logit(phi[i,t]) <- beta_phi1 + beta_phi2*elev[i] + rand[i] + trand[t] #if you wanted to simulate for a new time point, would be adding some noise to survival probability. You can say time effect is 0 for new time points to assume it's average.But makes more sense to just let JAGs do it, gives more uncertainty in following year but that makes sense. Makes sense for sites too.
  
  
  # Flowering:
  logit(kappa[i,t]) = beta_f[1] + beta_f[2] * flo[i,t] 
  
  # Relative abundance
  logit(lambda[i,t]) = beta_a[1] + beta_a[2] * dem[i,t]

  
  #  Pairwise 'source strength' calc and colonisation probability
  for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
  gammaDistPairs[i,n,t] <- 1-gamma0*exp(-gamma0*D.nb[i,n])*z[NB.mat[i,n],t] #This indexes the neighbour in question (NB.mat[i,n])
  }
  
  # Colonization prob
  gamma[i,t] = 1 - prod(1-gammaDistPairs[i,1:n.nb[i],t])
  
  
  ## Pr(z=1|z=1) = 1 - Pr(extinction)*Pr(not colonized) = 1 - ((1-phi)*(1-gamma))
  Ez[i,t+1] = (1 - z[i,t])*gamma[i,t] + (1-(1-phi[i,t])*(1-gamma[i,t]))*z[i,t] 
  #For year 2, colonization if it was not occupied and survival if it was alive the previous year 
  
  #new occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)
  p ~ dbeta(1, 1)
  
  #For random effects of year and time
  taut ~ dgamma(0.1,0.1)
  tau ~ dgamma(0.1,0.1)
  
  #For missing data in flowering and dem data
  mu_flo ~ dnorm(0, 0.001) #won't be negative! (maybe use a total vs. average <- poisson for the respoonse)
  tau_flo ~ dgamma(0.001,0.001) #so could use the total counts, use dgamma as conjugate to poisson for mean
  mu_dem ~ dgamma(0.001, 0.001)
  tau_dem ~ dgamma(0.001, 0.001)
  
  #conugate for mean paramater for normal dist is dnorm, for the precision its dgamma.
  #For predictors
  ##Survival
  beta_phi1 ~ dnorm(0, 1/1000)
  beta_phi2 ~ dnorm(0,1/1000)
  
  ##Colonization
  gamma0 ~ dbeta(1,1)
  
  ###Flowering
  beta_f[1] ~ dnorm(0,1/1000)
  beta_f[2] ~ dnorm(0, 1/1000)
  
  ###Relative abundance
  beta_a[1] ~ dnorm(0,1/1000)
  beta_a[2] ~ dnorm(0,1/1000)
  
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

Data_simple <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z, dem=dem, flo=flo, elev=elev[,1])

# 3) Specify a function to generate inital values for the parameters
#init.pa = occ
inits_fn = function() list(gamm0=1, psi1 = 0.1, mu_dem = 10, tau_dem= 10^-3, mu_flo = 10, tau_flo= 10^-3, beta_phi1=2, beta_phi2=2, alpha=2, z = z, p = 0.5)
#inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 2, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names = c('n_occ', 'p', 'gamma0', 'beta_phi1', 'beta_phi2', 'beta_a[1]', 'beta_a[2]', 'beta_f[1]', 'beta_f[2]', 'tau', 'taut')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)



#### MODEL 9 ####

#As before but includes a different distance decay

#Create list of neighbours within max distance (using 500 for now)
#Create a list of neighbours within max.dist for each patch
max.dist <- 200
NB.list <- apply(distmat,1, function(x) which(x <= max.dist))
# Number of neighbours per patch
n.nb <- sapply(NB.list,length)
summary(n.nb)
# 
# Here I create a matrix NB.mat with the patch numbers (indices) of all neighbours
# for each patch
# and a re-formated distance matrix D.nb with the distances to these neighbours
NB.mat <- matrix(NA, nrow = nrow(occ), ncol = max(n.nb))
D.nb <- NB.mat
for(i in 1:nrow(occ)){
  NB.mat[i,1:n.nb[i]] <- NB.list[[i]]
  D.nb[i,1:n.nb[i]] <- distmat[i,NB.list[[i]]]
}


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
  
  for(t in 1:(t.max-1)){
  trand[t] ~ dnorm(0, taut) } #if a time point is bad, it's bad for all sites so keep additive
  
  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)
  
  for(t in 1:(t.max-1)){
  
  #Model of missing data for dem and flo (in t loop because needs to be defined multiple times
  flo[i,t] ~ dnorm(mu_flo, tau_flo) #dpois if total, may likely need a negative binomial (search help for dnbin expects number of failures before successes, need to convert it). Or if you have enough observations, you could assume this is normal and the variances are estimated independently (vs. the poisson)
  dem[i,t] ~ dnorm(mu_dem, tau_dem) 
  
  #  Model of local survival (1-extinction) at site i
  logit(phi[i,t]) <- beta_phi1 + beta_phi2*elev[i] + trand[t] #if you wanted to simulate for a new time point, would be adding some noise to survival probability. You can say time effect is 0 for new time points to assume it's average.But makes more sense to just let JAGs do it, gives more uncertainty in following year but that makes sense. Makes sense for sites too.
  
  
  # Flowering:
  logit(kappa[i,t]) = beta_f[1] + beta_f[2] * flo[i,t] 
  
  # Relative abundance
  logit(lambda[i,t]) = beta_a[1] + beta_a[2] * dem[i,t]
  
  
  #  Pairwise 'source strength' calc and colonisation probability
  for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
  gammaDistPairs[i,n,t] <- 1-gamma0*exp(-D.nb[i,n]^2/(2*sigma^2))*z[NB.mat[i,n],t] 
  #This indexes the neighbour in question (NB.mat[i,n])
  #sigma is the distance decay scale parameter
  #gamma0 is baseline colonisation probability for coincident sites
  }
  
  # Colonization prob
  gamma[i,t] = 1 - prod(1-gammaDistPairs[i,1:n.nb[i],t]) #really this is site-level 'connectivity'
  #Could generate this at every site to look at connectivity at the landscape level, or sum to get a landscape scale connectivity measure
  
  
  ## Pr(z=1|z=1) = 1 - Pr(extinction)*Pr(not colonized) = 1 - ((1-phi)*(1-gamma))
  Ez[i,t+1] = gamma[i,t]*(1 - z[i,t]) + (1-(1-phi[i,t])*(1-gamma[i,t]))*z[i,t] 
  #P(occ) = P(colonized if not there) + 1-(P(extinction)*(P(not colonized and there, i.e. survived from last year)))
  #P(occ) = P(colonisation if not there) + P(survival if occupied already)
  #For year 2, colonization if it was not occupied and survival if it was alive the previous year 
  
  #new occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)
  p ~ dbeta(1, 1)
  
  #For random effects of year 
  taut ~ dgamma(0.1,0.1)
  
  #For missing data in flowering and dem data
  mu_flo ~ dnorm(0, 0.001) #won't be negative! (maybe use a total vs. average <- poisson for the respoonse)
  tau_flo ~ dgamma(0.001,0.001) #so could use the total counts, use dgamma as conjugate to poisson for mean
  mu_dem ~ dgamma(0.001, 0.001)
  tau_dem ~ dgamma(0.001, 0.001)
  
  #conugate for mean paramater for normal dist is dnorm, for the precision its dgamma.
  #For predictors
  ##Survival
  beta_phi1 ~ dnorm(0, 1/1000)
  beta_phi2 ~ dnorm(0,1/1000)
  
  ##Colonization
  gamma0 ~ dbeta(1,1)
  sigma ~ dgamma(1, 0.1)
  
  ###Flowering
  beta_f[1] ~ dnorm(0,1/1000)
  beta_f[2] ~ dnorm(0, 1/1000)
  
  ###Relative abundance
  beta_a[1] ~ dnorm(0,1/1000)
  beta_a[2] ~ dnorm(0,1/1000)
  
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

Data_simple <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z, dem=dem, flo=flo, elev=elev[,1])

# 3) Specify a function to generate inital values for the parameters
#init.pa = occ
inits_fn = function() list(sigma=10, gamma0=0.1, psi1 = 0.1, mu_dem = 1, tau_dem= 0.01, mu_flo = 1, tau_flo= 0.01, beta_phi1=1, beta_phi2=1, z = z, p = 0.9)
#inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 2, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names = c('n_occ', 'p', 'gamma0', 'sigma', 'beta_phi1', 'beta_phi2', 'beta_a[1]', 'beta_a[2]', 'beta_f[1]', 'beta_f[2]', 'taut')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)




#### MODEL 10 ####

#Model 9 but removes useless flower and abundance (where to put these??)

#Create list of neighbours within max distance (using 500 for now)
#Create a list of neighbours within max.dist for each patch
max.dist <- 200
NB.list <- apply(distmat,1, function(x) which(x <= max.dist))
# Number of neighbours per patch
n.nb <- sapply(NB.list,length)
summary(n.nb)
# 
# Here I create a matrix NB.mat with the patch numbers (indices) of all neighbours
# for each patch
# and a re-formated distance matrix D.nb with the distances to these neighbours
NB.mat <- matrix(NA, nrow = nrow(occ), ncol = max(n.nb))
D.nb <- NB.mat
for(i in 1:nrow(occ)){
  NB.mat[i,1:n.nb[i]] <- NB.list[[i]]
  D.nb[i,1:n.nb[i]] <- distmat[i,NB.list[[i]]]
}


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
  
  for(t in 1:(t.max-1)){
  trand[t] ~ dnorm(0, taut) } #if a time point is bad, it's bad for all sites so keep additive
  
  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)
  
  for(t in 1:(t.max-1)){
  
  #Model of missing data for dem and flo (in t loop because needs to be defined multiple times
  flo[i,t] ~ dnorm(mu_flo, tau_flo) #dpois if total counts (norm if an average)
  dem[i,t] ~ dnorm(mu_dem, tau_dem) 
  
  #  Model of local survival (1-extinction) at site i
  logit(phi[i,t]) <- beta_phi1 + beta_phi2*elev[i] + trand[t] #if you wanted to simulate for a new time point, would be adding some noise to survival probability. You can say time effect is 0 for new time points to assume it's average.But makes more sense to just let JAGs do it, gives more uncertainty in following year but that makes sense. Makes sense for sites too.
  
  #  Pairwise 'source strength' calc and colonisation probability
  for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
  gammaDistPairs[i,n,t] <- 1-gamma0*exp(-(D.nb[i,n]^2)/(2*sigma^2))*z[NB.mat[i,n],t] 
  #This indexes the neighbour in question (NB.mat[i,n])
  #sigma is the distance decay scale parameter
  #gamma0 is baseline colonisation probability for coincident sites
  }
  
  # Colonization prob
  gamma[i,t] = 1 - prod(1-gammaDistPairs[i,1:n.nb[i],t]) #really this is site-level 'connectivity'
  #Could generate this at every site to look at connectivity at the landscape level, or sum to get a landscape scale connectivity measure
  
  
  ## Pr(z=1|z=1) = 1 - Pr(extinction)*Pr(not colonized) = 1 - ((1-phi)*(1-gamma))
  Ez[i,t+1] = gamma[i,t]*(1 - z[i,t]) + (1-(1-phi[i,t])*(1-gamma[i,t]))*z[i,t] 
  #P(occ) = P(colonized if not there) + 1-(P(extinction)*(P(not colonized and there, i.e. survived from last year)))
  #P(occ) = P(colonisation if not there) + P(survival if occupied already)
  #For year 2, colonization if it was not occupied and survival if it was alive the previous year 
  
  #new occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)
  p ~ dbeta(1, 1)
  
  #For random effects of year 
  taut ~ dgamma(0.1,0.1)
  
  #conugate for mean paramater for normal dist is dnorm, for the precision its dgamma.
  #For predictors
  ##Survival
  beta_phi1 ~ dnorm(0, 1/1000)
  beta_phi2 ~ dnorm(0,1/1000)
  
  ##Colonization
  gamma0 ~ dbeta(1,1)
  sigma ~ dgamma(1, 0.1)
  
  ###Derived quantities block
  psi[1]=psi1
  n_occ[1]=sum(z[1:n.sites,1])

  for(t in 1:(t.max-1)) {
  n_occ[t+1]=sum(z[1:n.sites,t+1])
  }
  
  # for (t in 1:(t.max-1)) {
  # n_gamma[t+1] = sum(gamma[1:n.sites,t+1])
  # }

  
  
  }
  ",fill = TRUE, file="occ1.txt")

sink()
#you could ask just for the mean, is in rjags! Look at means only, not posterior (example in WAIC folder in Day 3)
#Put in truncated dispersal kernel, estimate range limit


### 2) Set up a list that contains all the necessary data

Data_simple <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z, elev=elev[,1])

# 3) Specify a function to generate inital values for the parameters
#init.pa = occ
inits_fn = function() list(sigma=2, gamma0=0.1, psi1 = 0.1, beta_phi1=1, beta_phi2=1, z = z, p = 0.9)
#inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 3, n.adapt= 2000)
# Specify parameters for which posterior samples are saved
para.names = c('n_occ', 'p', 'gamma0', 'sigma', 'beta_phi1', 'beta_phi2', 'gamma', 'taut')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 2000)

par_vals = summary(Samples)$statistics
summary(Samples)
#plot(Samples,ask=T)


#### MODEL 11 ####

#Model 9, using negative exponential and weighting by patch abundance and flowering

#Create list of neighbours within max distance (using 500 for now)
#Create a list of neighbours within max.dist for each patch
max.dist <- 200
NB.list <- apply(distmat,1, function(x) which(x <= max.dist))
# Number of neighbours per patch
n.nb <- sapply(NB.list,length)
summary(n.nb)
# 
# Here I create a matrix NB.mat with the patch numbers (indices) of all neighbours
# for each patch
# and a re-formated distance matrix D.nb with the distances to these neighbours
NB.mat <- matrix(NA, nrow = nrow(occ), ncol = max(n.nb))
D.nb <- NB.mat
for(i in 1:nrow(occ)){
  NB.mat[i,1:n.nb[i]] <- NB.list[[i]]
  D.nb[i,1:n.nb[i]] <- distmat[i,NB.list[[i]]]
}


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
  
  for(t in 1:(t.max-1)){
  trand[t] ~ dnorm(0, taut) } #if a time point is bad, it's bad for all sites so keep additive
  
  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)
  
  for(t in 1:(t.max-1)){
  
  #Model of missing data for dem and flo (in t loop because needs to be defined multiple times
  flo[i,t] ~ dnorm(mu_flo, tau_flo) #dpois if total counts (norm for average)
  dem[i,t] ~ dnorm(mu_dem, tau_dem) 
  
  # Model of relative abundance
  logit(lambda[i,t]) = beta_a[1] + beta_a[2]*dem[i,t] 
  
  # Model of flowering
	logit(kappa[i,t]) = beta_f[1] + beta_f[2]*flo[i,t] 

  #  Model of local survival (1-extinction) at site i
  logit(phi[i,t]) <- beta_phi[1] + beta_phi[2]*gdd[i] + beta_phi[3]*fdd[i] + trand[t] #if you wanted to simulate for a new time point, would be adding some noise to survival probability. You can say time effect is 0 for new time points to assume it's average.But makes more sense to just let JAGs do it, gives more uncertainty in following year but that makes sense. Makes sense for sites too.
  
  #  Pairwise 'source strength' calc and colonisation probability
  for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
  gammaDistPairs[i,n,t] <- gamma0 * exp(-gamma0*D.nb[i,n]) * kappa[NB.mat[i,n],t] * lambda[NB.mat[i,n],t]
  }
  
  # Colonization prob
  gamma[i,t] = 1 - prod(1-gammaDistPairs[i,1:n.nb[i],t]) #really this is site-level 'connectivity'
  #Could generate this at every site to look at connectivity at the landscape level, or sum to get a landscape scale connectivity measure
  
  
  ## Pr(z=1|z=1) = 1 - Pr(extinction)*Pr(not colonized) = 1 - ((1-phi)*(1-gamma))
  Ez[i,t+1] = gamma[i,t]*(1 - z[i,t]) + (1-(1-phi[i,t])*(1-gamma[i,t]))*z[i,t] 
  #P(occ) = P(colonized if not there) + 1-(P(extinction)*(P(not colonized and there, i.e. survived from last year)))
  #P(occ) = P(colonisation if not there) + P(survival if occupied already)
  #For year 2, colonization if it was not occupied and survival if it was alive the previous year 
  
  #new occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)
  p ~ dbeta(1, 1)
  
  #For random effects of year 
  taut ~ dgamma(0.1,0.1)

  #For missing data in flowering and dem data
  mu_flo ~ dnorm(0, 0.001) #won't be negative! (maybe use a total vs. average <- poisson for the respoonse)
  tau_flo ~ dgamma(0.001,0.001) #so could use the total counts, use dgamma as conjugate to poisson for mean
  mu_dem ~ dnorm(0, 0.001)
  tau_dem ~ dgamma(0.001, 0.001)
  
  #conjugate for mean paramater for normal dist is dnorm, for the precision its dgamma.
  #For predictors
  ##Survival
  beta_phi[1] ~ dnorm(0, 1/1000)
  beta_phi[2] ~ dnorm(0,1/1000)
  beta_phi[3] ~ dnorm(0,1/1000)

  ##Flowering
  beta_f[1] ~ dnorm(0,1/1000)
  beta_f[2] ~ dnorm(0, 1/1000)
  
  ##Relative abundance
  beta_a[1] ~ dnorm(0,1/1000)
  beta_a[2] ~ dnorm(0,1/1000)
  
  ##Colonization
  gamma0 ~ dbeta(1,1)
  
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

Data_simple <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z, dem=dem, flo=flo, fdd=fdd[,1], gdd=gdd[,1])


# 3) Specify a function to generate inital values for the parameters
inits_fn = function() list(gamma0=0.1, psi1 = 0.1, mu_dem = 1, tau_dem= 0.01, mu_flo = 1, tau_flo= 0.01, beta_phi=runif(3,-3,3), beta_a=runif(2, -3,3), beta_f=runif(2,-3,3), z = z, p = 0.9)

load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 1, n.adapt= 1000)

# Specify parameters for which posterior samples are saved
para.names = c('n_occ', 'p', 'gamma0', 'beta_phi[1]', 'beta_phi[2]', 'beta_phi[3]', 'beta_a[1]', 'beta_a[2]', 'beta_f[1]', 'beta_f[2]', 'taut')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)






#### MODEL 12 ####

#Model 11, using negative exponential and weighting by patch abundance and flowering (plus climate variables)

#Create list of neighbours within max distance (using 500 for now)
#Create a list of neighbours within max.dist for each patch
max.dist <- 200
NB.list <- apply(distmat,1, function(x) which(x <= max.dist))
# Number of neighbours per patch
n.nb <- sapply(NB.list,length)
summary(n.nb)
# 
# Here I create a matrix NB.mat with the patch numbers (indices) of all neighbours
# for each patch
# and a re-formated distance matrix D.nb with the distances to these neighbours
NB.mat <- matrix(NA, nrow = nrow(occ), ncol = max(n.nb))
D.nb <- NB.mat
for(i in 1:nrow(occ)){
  NB.mat[i,1:n.nb[i]] <- NB.list[[i]]
  D.nb[i,1:n.nb[i]] <- distmat[i,NB.list[[i]]]
}


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
  
  #Random year effect (for survival)
  for(t in 1:(t.max-1)){
  trand[t] ~ dnorm(0, taut) }
  
  # State (process) model
  # 1st year (2008) 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)
  
  # Years 2009-2012
  for(t in 1:(t.max-1)){
  
  #Modelling missing data for dem and flo
  flo[i,t] ~ dnorm(mu_flo, tau_flo) #dpois if total counts (norm for average)
  dem[i,t] ~ dnorm(mu_dem, tau_dem) 
  
  # Model of relative abundance
  logit(lambda[i,t]) = beta_a[1] + beta_a[2]*dem[i,t] + beta_a[3]*gdd[i]
  
  # Model of flowering
  logit(kappa[i,t]) = beta_f[1] + beta_f[2]*flo[i,t] 
  
  #  Pairwise 'source strength' calc and colonisation probability
  for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
  gammaDistPairs[i,n,t] <- gamma0 * exp(-gamma0*D.nb[i,n]) * kappa[NB.mat[i,n],t] * lambda[NB.mat[i,n],t] * z[NB.mat[i,n],t]
  }
  
  # Model of colonization probability 
  gamma[i,t] = 1 - prod(1-gammaDistPairs[i,1:n.nb[i],t]) #really this is site-level 'connectivity'
  
  #  Model of local survival probability (1-extinction) 
  logit(phi[i,t]) <- beta_phi[1] + beta_phi[2]*gdd[i] + beta_phi[3]*fdd[i] + trand[t] 

  # Generating occupancy probability
  Ez[i,t+1] = gamma[i,t]*(1 - z[i,t]) + (1-(1-phi[i,t])*(1-gamma[i,t]))*z[i,t] 
  #P(occ) = P(colonized if not there) + 1-(P(extinction)*(P(not colonized and there, i.e. survived from last year)))
  
  #True occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  }
  }
  
  
  # Prior distributions
  ##For year 1
  psi1 ~ dbeta(1,1)
  p ~ dbeta(1, 1)
  
  #For random effects of year 
  taut ~ dgamma(0.1,0.1)
  
  #For missing data in flowering and dem data
  mu_flo ~ dnorm(0, 0.001) #won't be negative! (maybe use a total vs. average <- poisson for the respoonse)
  tau_flo ~ dgamma(0.001,0.001) #so could use the total counts, use dgamma as conjugate to poisson for mean
  mu_dem ~ dnorm(0, 0.001)
  tau_dem ~ dgamma(0.001, 0.001)
  
  #conjugate for mean paramater for normal dist is dnorm, for the precision its dgamma.
  #For predictors
  ##Survival
  beta_phi[1] ~ dnorm(0, 1/1000)
  beta_phi[2] ~ dnorm(0,1/1000)
  beta_phi[3] ~ dnorm(0,1/1000)
  
  ##Flowering
  beta_f[1] ~ dnorm(0,1/1000)
  beta_f[2] ~ dnorm(0, 1/1000)
  
  ##Relative abundance
  beta_a[1] ~ dnorm(0,1/1000)
  beta_a[2] ~ dnorm(0,1/1000)
  beta_a[3] ~ dnorm(0,1/1000)
  
  ##Colonization
  gamma0 ~ dbeta(1,1)
  
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

Data_simple <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z, dem=dem, flo=flo, fdd=fdd[,1], gdd=gdd[,1])


# 3) Specify a function to generate inital values for the parameters
inits_fn = function() list(gamma0=0.1, psi1 = 0.1, mu_dem = 1, tau_dem= 0.01, mu_flo = 1, tau_flo= 0.01, beta_phi=runif(3,-3,3), beta_a=runif(3, -3,3), beta_f=runif(2,-3,3), z = z, p = 0.9)

load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 1, n.adapt= 1000)

# Specify parameters for which posterior samples are saved
para.names = c('n_occ', 'p', 'gamma0', 'beta_phi[1]', 'beta_phi[2]', 'beta_phi[3]', 'beta_a[1]', 'beta_a[2]', 'beta_a[2]', 'beta_f[1]', 'beta_f[2]', 'taut')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)






#### MODEL 13 ####

#Model 11, using negative exponential and weighting by patch abundance and flowering (plus climate variables)

# Load data
load('jags_data.RData')

# Create list of neighbours within max distance (using 500 for now)
max.dist <- 500
NB.list <- apply(jags_data$distmat,1, function(x) which(x <= max.dist))


# Number of neighbours per patch
n.nb <- sapply(NB.list,length)
summary(n.nb)

# A re-formated distance matrix D.nb with the distances to these neighbours
NB.mat <- matrix(NA, nrow = nrow(jags_data$occ), ncol = max(n.nb))
D.nb <- NB.mat
for(i in 1:nrow(jags_data$occ)){
  NB.mat[i,1:n.nb[i]] <- NB.list[[i]]
  D.nb[i,1:n.nb[i]] <- distmat[i,NB.list[[i]]]
}

# Scaling predictors
#library(scales) #this is for rescale, which changes it to between 0 and 1, but I already create
#a plogit transformation below
#NEED TO STANDARDIZE SCALE ACROSS MATRIX!!!!!! ALSO NEED TO STATE PREDICTORS NEEDED TO STAY ABOVE 0.
dem <- as.matrix(scale(jags_data$dem))
flo <- as.matrix(scale(jags_data$flo))
pat <- as.matrix(scale(jags_data$pat))
elev <- c(scale(jags_data$elev[,1]))


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
  # 1st year (2008) 
  for(i in 1:n.sites){
  z[i,1] ~ dbern(psi1)
  
  # Years 2009-2012
  for(t in 1:(t.max-1)){
  
  #Modelling missing data for demographic data
  dem[i,t] ~ dnorm(mu_dem[i,t], tau_dem) 
  flo[i,t] ~ dnorm(mu_flo[i,t], tau_flo) 
  pat[i,t] ~ dnorm(mu_pat[i,t], tau_pat) 

  #Models of elevation effects on dem data
  mu_dem[i,t] <- beta_dem[1] + beta_dem[2]*elev[i]
  mu_flo[i,t] <- beta_flo[1] + beta_flo[2]*elev[i]
  mu_pat[i,t] <- beta_pat[1] + beta_pat[2]*elev[i]

  # Model of relative abundance
  logit(totseed[i,t]) <- mu_dem[i,t]*mu_flo[i,t]*mu_pat[i,t]
  
  #  Pairwise 'source strength' calc and colonisation probability
  for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
  gammaDistPairs[i,n,t] <- gamma0 * exp(-gamma0*D.nb[i,n]) * totseed[NB.mat[i,n],t] * z[NB.mat[i,n],t] 
  }

  # Model of colonization probability 
  gamma[i,t] = 1 - prod(1-gammaDistPairs[i,1:n.nb[i],t]) #really this is site-level 'connectivity'


  #  Model of local survival probability (1-extinction) 
  logit(phi[i,t]) <- beta_phi[1] + beta_phi[2] * elev[i] 
  
  # Generating occupancy probability
  Ez[i,t+1] = gamma[i,t]*(1 - z[i,t]) + (1-(1-phi[i,t])*(1-gamma[i,t]))*z[i,t] + 0.001
  #P(occ) = P(colonized if not there) + 1-(P(extinction)*(P(not colonized and there, i.e. survived from last year)))
  
  #True occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  
  }
  }
  
  
  # Prior distributions
  
  ##For year 1
  psi1 ~ dbeta(1,1)
  
  #Detection
  p ~ dbeta(1, 1)
  
  
  #For missing data in flowering and dem data
  #mu_dem ~ dnorm(0, 0.001) 
  tau_dem ~ dgamma(0.001,0.001) 
  #mu_flo ~ dnorm(0, 0.001) 
  tau_flo ~ dgamma(0.001,0.001) 
  #mu_pat ~ dnorm(0, 0.001) 
  tau_pat ~ dgamma(0.001,0.001) 
  
  #Demography
  beta_dem[1] ~ dnorm(0, 1/1000)
  beta_dem[2] ~ dnorm(0, 1/1000)
  beta_flo[1] ~ dnorm(0, 1/1000)
  beta_flo[2] ~ dnorm(0, 1/1000)
  beta_pat[1] ~ dnorm(0, 1/1000)
  beta_pat[2] ~ dnorm(0, 1/1000)

  #Survival
  beta_phi[1] ~ dnorm(0, 1/1000)
  beta_phi[2] ~ dnorm(0, 1/1000)
  
  #Colonization
  gamma0 ~ dbeta(1,1)
  
  ###Derived quantities block
  #First year
  psi[1]=psi1
  n_occ[1]=sum(z[1:n.sites,1])
  
  #All subsequent years
  for(t in 1:(t.max-1)) {
  n_occ[t+1]=sum(z[1:n.sites,t+1])
  }
  
  
  
  }
  ",fill = TRUE, file="occ1.txt")

sink()

### 2) Set up a list that contains all the necessary data

Data_simple <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(jags_data$occ), t.max = ncol(jags_data$occ), y = jags_data$occ, z=jags_data$z, dem=dem, flo=flo, pat=pat, elev=elev)

init.z <- as.matrix(jags_data$z)

init.z[is.na(init.z)] <- 0

# 3) Specify a function to generate inital values for the parameters
inits_fn = function() list(gamma0=0.1, psi1 = 0.1, tau_dem= 1, tau_flo= 1, tau_pat = 1, beta_dem=runif(2,-3,3), beta_flo=runif(2,-3,3), beta_pat=runif(2,-3,3), beta_phi=runif(2,-3,3), p = 0.5)

jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 1, n.adapt= 10, inits = inits_fn)

# Specify parameters for which posterior samples are saved
para.names = c('n_occ', 'p', 'gamma0', 
               'beta_dem[1]', 'beta_dem[2]', 
               'beta_flo[1]', 'beta_flo[2]', 
               'beta_pat[1]', 'beta_pat[2]',
               'beta_phi[1]', 'beta_phi[2]') 

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1500)

save(Samples, file='Sample_outputv2.RData')

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)
