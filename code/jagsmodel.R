#### Full model 14.01.2020

# Load data
load('jags_data.RData')

# Create list of neighbours within max distance (using 500 for now)
max.dist <- 200
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
library(scales)
dem <- rescale(jags_data$dem)
flo <- rescale(jags_data$flo)
pat <- rescale(jags_data$pat)
elev <- rescale(jags_data$elev[,1])


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
  
  #Modelling missing data for demographic data
  dem[i,t] ~ dnorm(mu_dem, tau_dem) 
  flo[i,t] ~ dnorm(mu_flo, tau_flo) 
  pat[i,t] ~ dnorm(mu_pat, tau_pat) 

  # Model of relative abundance
  totseed[i,t] <- dem[i,t]*flo[i,t]*pat[i,t]
  
  #  Pairwise 'source strength' calc and colonisation probability
  for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
  gammaDistPairs[i,n,t] <- gamma0 * exp(-gamma0*D.nb[i,n]) * totseed[NB.mat[i,n],t] * z[NB.mat[i,n],t] 
  }
  
  # Model of colonization probability 
  gamma[i,t] = 1 - prod(1-gammaDistPairs[i,1:n.nb[i],t]) #really this is site-level 'connectivity'

  #  Model of local survival probability (1-extinction) + random year effect 
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

  #Detection
  p ~ dbeta(1, 1)
  
  #For random effects of year 
  taut ~ dgamma(0.001,0.001)
  
  #For missing data in flowering and dem data
  mu_dem ~ dnorm(0, 0.001) 
  tau_dem ~ dgamma(0.001,0.001) 
  mu_flo ~ dnorm(0, 0.001) 
  tau_flo ~ dgamma(0.001,0.001) 
  mu_pat ~ dnorm(0, 0.001) 
  tau_pat ~ dgamma(0.001,0.001) 
  
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
#you could ask just for the mean, is in rjags! Look at means only, not posterior (example in WAIC folder in Day 3)
#Put in truncated dispersal kernel, estimate range limit


### 2) Set up a list that contains all the necessary data

Data_simple <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(jags_data$occ), t.max = ncol(jags_data$occ), y = jags_data$occ, z=jags_data$z, dem=dem, flo=flo, pat=pat, elev=elev)


# 3) Specify a function to generate inital values for the parameters
inits_fn = function() list(gamma0=0.1, psi1 = 0.1, mu_dem = 1, tau_dem= 1, mu_flo = 1, tau_flo= 1, mu_pat = 1, tau_pat = 1, beta_phi=runif(2,-3,3), z = z, p = 0.9)

jagsModel = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 3, n.adapt= 50)

# Specify parameters for which posterior samples are saved
para.names = c('n_occ', 'p', 'gamma0', 'beta_phi[1]', 'beta_phi[2]', 'taut') 

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)


