#### Solidago dynamic occupancy model ####
# C. Chisholm, 04.2021
# This model estimates occupancy of Solidago canadensis over five years of data (2008-2012, except no data is available for 2009) in the Rhine valley (eastern Swiss alps).
# Occupancy data is available along a transect across the valley bottom (following cycleways and roads).
# Also availale are data on population density and flowering number for a subset of patches along this transect in years 2010,2011 and 2012.
# The following model estimates colonization probability (gamma) based on the distance scaled input of seed number (totseed, obtained from pop density N and individual flowering number flo).
# Gamma is also a function of elevation (through effects on seed production as well as a direct effect on gamma, assumed due to an effect of climate on establishment)
# Survival (phi) is a function of elevation.

library(runjags)

### 1) Load data and transform 

load('jags_data.RData')

# This is to create a subset of 500 patches, remove if running full dataset
# jags_data$occ <- jags_data$occ[1:500,]
# jags_data$z <- jags_data$z[1:500,]
# jags_data$N <- jags_data$N[1:500,]
# jags_data$flo <- jags_data$flo[1:500,]
# jags_data$elev <- jags_data$elev[1:500]
# jags_data$distmat <- jags_data$distmat[1:500,1:500]

# Create list of neighbours within max distance (using 1000 mfor now)
max.dist <- 1000
NB.list <- apply(jags_data$distmat,1, function(x) which(x <= max.dist))

# Number of neighbours per patch
n.nb <- sapply(NB.list,length)
summary(n.nb)

# A re-formated distance matrix D.nb with the distances to these neighbours
NB.mat <- matrix(NA, nrow = nrow(jags_data$occ), ncol = max(n.nb))
D.nb <- NB.mat
for(i in 1:nrow(jags_data$occ)){
  NB.mat[i,1:n.nb[i]] <- NB.list[[i]]
  D.nb[i,1:n.nb[i]] <- jags_data$distmat[i,NB.list[[i]]]
}

# Scaling predictors betwwen 0 and 1
library(scales)
elev <- scale(jags_data$elev)[,]

#Remove 4 cells with 0 values in N and flo
jags_data$N[jags_data$N==0]<-NA
jags_data$flo[jags_data$flo==0]<-NA

jags_data$logN <- log(jags_data$N)
jags_data$logflo <- log(jags_data$flo)

### 2) Occupancy model

sink("occ1.txt") 

cat(
  "
  model{
  
  ## Observation model
  for(t in 1:t.max){
  for(i in 1:n.sites){
  muY[i,t] <- z[i,t]*p
  y[i,t] ~ dbern(muY[i,t])
  }
  }
  
  #Random year effect (for pop parameteres)
  #for(t in 1:(t.max-1)){
  #trand[t] ~ dnorm(0, taut) 
  #}
  
  ## State (process) model
  # 1st year (2008) 
  for(i in 1:n.sites){
  Ez[i,1] <- psi1
  z[i,1] ~ dbern(Ez[i,1])
  }
  
  # Years 2009-2012
  for (t in 1:(t.max-1)){
  for(i in 1:n.sites){
  
  # Modelling missing data for demographic data (data scaled!--> zeros so use lognormal to constrain to positive values)
  N[i,t] ~ dnorm(mu_N[i,t], tau_N)
  flo[i,t] ~ dnorm(mu_flo[i,t], tau_flo)
  
  mu_N[i,t] <- mu_n[i,t]*z[i,t]
  mu_flo[i,t] <- mu_f[i,t]*z[i,t]
  log(mu_n[i,t]) <- beta_N[1] + beta_N[2]*elev[i] 
  log(mu_f[i,t]) <- beta_flo[1] + beta_flo[2]*elev[i] 
  
  # Creating neighboring patch weights [0-1] by total seed output per patch
  logit(totseed[i,t]) <- N[i,t]*flo[i,t]
  
  # Neighbourhood connectivity to focal patch i (weighted by total seed production of neighbourhood)
  for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
  gammaDistPairs[i,n,t] <- totseed[NB.mat[i,n],t] * exp(-alpha*D.nb[i,n]) * z[NB.mat[i,n],t]
  }
  
  # Source strength calculation per patch
  S[i,t] = prod(1-gammaDistPairs[i,1:n.nb[i],t]) #site-level 'connectivity'
  
  # Model of colonization probability 
  gamma[i,t] = 1 - exp(-rho[i,t]*S[i,t])
  
  logit(rho[i,t]) <- beta_rho[1] + beta_rho[2]*elev[i] 
  
  # Model of local survival probability (1-extinction) 
  logit(phi[i,t]) <- beta_phi[1] + beta_phi[2]*elev[i] 
  
  # Generating occupancy probability
  Ez[i,t+1] = gamma[i,t]*(1 - z[i,t]) + phi[i,t]*z[i,t]
  #P(occ) = P(colonized if not there) + 1-(P(extinction)*(P(not colonized and there, i.e. survived from last year)))
  
  # True occupancy probability  
  z[i,t+1] ~ dbern(Ez[i,t+1])
  
  }
  }
  
  
  ### Prior distributions
  
  #For year 1
  psi1 ~ dbeta(1,1)
  
  # For detection probability
  p ~ dbeta(480, 271) #based on 2020 data
  
  #For random effects of year 
  #taut ~ dgamma(0.001,0.001)
  
  #Demography
  beta_N[1] ~ dnorm(0, 0.37)
  beta_N[2] ~ dnorm(0, 0.37)
  tau_N ~ dexp(0.8)  #dgamma(2,1) works well
  
  beta_flo[1] ~ dnorm(0, 0.37)
  beta_flo[2] ~ dnorm(0, 0.37)
  tau_flo ~ dexp(0.8)  #was dgamma(10,10) (worked great), 0.1 was ok

  #Colonisation
  beta_rho[1] ~ dnorm(0, 0.37)
  beta_rho[2] ~ dnorm(0, 0.37)
  
  #Survival
  beta_phi[1] ~ dnorm(0, 0.37) 
  beta_phi[2] ~ dnorm(0, 0.37)
  
  #Colonization
  alpha ~ dexp(1) #dispersal scale parameter
  
  ###Derived quantities block
  #First year
  psi[1] <- psi1
  n_occ[1]<- sum(z[1:n.sites,1])
  
  #All subsequent years
  for(t in 1:(t.max-1)) {
  n_occ[t+1] <- sum(z[1:n.sites,t+1])
  }
  
  }
  ",
  fill = TRUE, file="occ1.txt")

sink()

### 3) Set up a list that contains all the necessary data
dat <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(jags_data$occ), t.max = ncol(jags_data$occ), y = jags_data$occ, N=jags_data$logN, flo=jags_data$logflo, elev=elev) 

### 4) Create initial values
# Specify a function to generate initial values for the parameters
#inits_fn = function() list(z=matrix(data=rep(1, 500*5), nrow=500, ncol=5), alpha=rexp(1,1), psi1 = runif(1,0.01,1), tau_flo= rgamma(1, 0.1, 0.1), tau_N= rgamma(1, 0.1, 0.1), beta_N=rnorm(2, 0, 1.6), beta_flo=rnorm(2, 0, 1.6), beta_phi=rnorm(2,0,1.6), beta_rho=rnorm(2, 0, 1.6)) 
inits_fn = function() list(z=matrix(data=rep(1, 2741), nrow=2741, ncol=5), alpha=rexp(1,1), psi1 = runif(1,0.01,1), tau_flo= rgamma(1, 0.1, 0.1), tau_N= rgamma(1, 0.1, 0.1), beta_N=rnorm(2, 0, 1.6), beta_flo=rnorm(2, 0, 1.6), beta_phi=rnorm(2,0,1.6), beta_rho=rnorm(2, 0, 1.6)) 

### 5) Specify parameters for which posterior samples are saved
para.names = c('psi1', 'n_occ', 'alpha', 'tau_flo', 'tau_N',
               'beta_rho', 'beta_phi', 'beta_N', 'beta_flo') 

### 6) Run model
parsamples <- run.jags('occ1.txt', data=dat, inits = inits_fn, n.chains = 3, monitor=para.names,
                       adapt = 1000, burnin = 10000, thin=30, sample = 10000, method='parallel')


### 7) Model diagnostics
load('modeloutput062021.RData')
summary(parsamples)
plot(parsamples)
gelman.diag(parsamples$mcmc)

