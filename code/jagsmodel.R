#### Full model 14.01.2020

#Model includes population predictors and climate
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

library(scales)
totseed <- rescale(dem*pa*flo)
elev <- rescale(elev
                )
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
  totseed[i,t] ~ dbeta(mu_seed, tau_seed) #dpois if total counts (norm for average)
  
  # Model of relative abundance
  totseed[i,t] = beta_a[1] + beta_a[2] * elev[i] 
  
  #  Pairwise 'source strength' calc and colonisation probability
  for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
  gammaDistPairs[i,n,t] <- gamma0 * exp(-gamma0*D.nb[i,n]) * totseed[NB.mat[i,n],t] * z[NB.mat[i,n],t] 
  }
  
  # Model of colonization probability 
  gamma[i,t] = 1 - prod(1-gammaDistPairs[i,1:n.nb[i],t]) #really this is site-level 'connectivity'

  #  Model of local survival probability (1-extinction) 
  logit(phi[i,t]) <- beta_phi[1] + beta_phi[2]*elev[i] + trand[t] 
  
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
  taut ~ dgamma(0.001,0.001)
  
  #For missing data in flowering and dem data
  mu_seed ~ dnorm(0, 0.001) #shape1
  tau_seed ~ dnorm(0.001,0.001) #shape2
  
  #conjugate for mean paramater for normal dist is dnorm, for the precision its dgamma.
  #For predictors
  ##Survival
  beta_phi[1] ~ dnorm(0, 1/1000)
  beta_phi[2] ~ dnorm(0,1/1000)
  
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

Data_simple <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(occ), t.max = ncol(occ), y = occ, z=z, totseed=totseed, elev=elev[,1])


# 3) Specify a function to generate inital values for the parameters
inits_fn = function() list(gamma0=0.1, psi1 = 0.1, mu_seed = 1, tau_seed= 1, beta_phi=runif(2,-3,3), beta_a=runif(2, -3,3),z = z, p = 0.9)

load.module('glm')
jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 1, n.adapt= 1000)

# Specify parameters for which posterior samples are saved
para.names = c('n_occ', 'p', 'gamma0', 'beta_phi[1]', 'beta_phi[2]', 'beta_a[1]', 'beta_a[2]', 'beta_f[1]', 'beta_f[2]', 'taut')  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)


