#### Full model 30.03.2020

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
  D.nb[i,1:n.nb[i]] <- jags_data$distmat[i,NB.list[[i]]]
}

# Scaling predictors betwwen 0 and 1
library(scales)
elev <- scale(jags_data$elev)[,]
#N <- as.matrix(scale(jags_data$N))[,]
#flo <- as.matrix(scale(jags_data$flo))[,]
#elev <- (jags_data$elev-mean(jags_data$elev))/sd(jags_data$elev)
#anthro <- (jags_data$anthro-mean(jags_data$anthro))/sd(jags_data$anthro)

#Create NA columns to simulate data to 2020 (+8)
# jags_data$occ <- as.matrix(cbind(jags_data$occ, matrix(NA, 2741, 8)))
# jags_data$z <- as.matrix(cbind(jags_data$z, matrix(NA, 2741, 8)))
# jags_data$N <- as.matrix(cbind(jags_data$N, matrix(NA, 2741, 8)))
# jags_data$flo <- as.matrix(cbind(jags_data$flo, matrix(NA, 2741, 8)))
# jags_data$nyear <- 13

#load.module('glm')

sink("occ1.txt") 

cat(
  "
model{
  
  ## Observation model
  for(t in 1:t.max){
    for(i in 1:n.sites){
      muY[i,t] <- z[i,t]
      y[i,t] ~ dbern(muY[i,t])
    }
  }
  
  #Random year effect (for survival)
  for(t in 1:(t.max-1)){
    trand[t] ~ dnorm(0, taut) 
  }
  
  ## State (process) model
  # 1st year (2008) 
  for(i in 1:n.sites){
    Ez[i,1] <- psi1
    z[i,1] ~ dbern(Ez[i,1])
  }

  
  # Years 2009-2012
  for (t in 1:(t.max-1)){
    for(i in 1:n.sites){
  
  # Modelling missing data for demographic data (data scaled!)
  N[i,t] ~ dnorm(mu_N[i,t], tau_N) 
  flo[i,t] ~ dnorm(mu_flo[i,t], tau_flo) 

  # Models of elevation effects on dem data
  mu_N[i,t] <- beta_N[1] + beta_N[2]*elev[i] 
  mu_flo[i,t] <- beta_flo[1] + beta_flo[2]*elev[i] 
  
  # Model of relative abundance
  logit(totseed[i,t]) <- N[i,t]*flo[i,t]
  
  # Neighbourhood connectivity to focal patch i (weighted by total seed production of neighbourhood)
    for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
      gammaDistPairs[i,n,t] <- totseed[NB.mat[i,n],t] * exp(-alpha*D.nb[i,n]) * z[NB.mat[i,n],t]
    }
  
  # Source strength calc
  S[i,t] = prod(1-gammaDistPairs[i,1:n.nb[i],t]) #really this is site-level 'connectivity'

  # Model of colonization probability 
  gamma.z[i,t] = 1 - exp(-rho[i,t]*S[i,t])
  gamma[i,t]<-max(0.0001,min(0.9999, gamma.z[i,t]))

  logit(rho[i,t]) <- beta_rho[1] + beta_rho[2]*elev[i] + trand[t] 

  # Model of local survival probability (1-extinction) + random year effect 
  logit(phi.z[i,t]) <- beta_phi[1] + beta_phi[2]*elev[i] + trand[t] 
  phi[i,t]<-max(0.0001,min(0.9999, phi.z[i,t]))
  
  # Generating occupancy probability
  Ez[i,t+1] = gamma[i,t]*(1 - z[i,t]) + phi[i,t]*z[i,t] + 0.001
  #P(occ) = P(colonized if not there) + 1-(P(extinction)*(P(not colonized and there, i.e. survived from last year)))
  
  # True occupancy probability  
  z[i,t+1] ~ dbin(Ez[i,t+1],1)

    }
  }
  
  
  ### Prior distributions
  
  #For year 1
  psi1 ~ dbeta(1,1)

  # for (i in 1:n.sites) {
  #   for (t in 1:(t.max-1)) {
  #     gamma[i,t] ~ dunif(0,1)
  #     phi[i,t] ~ dunif(0,1)
  #   }
  # }
  
  #For random effects of year 
  taut ~ dgamma(0.001,0.001)
  
  #For missing data in flowering and dem data
  tau_N ~ dgamma(0.001,0.001) 
  tau_flo ~ dgamma(0.001,0.001) 
  
  #Demography
  beta_N[1] ~ dnorm(0, 1/1000)
  beta_N[2] ~ dnorm(0, 1/1000)
  beta_flo[1] ~ dnorm(0, 1/1000)
  beta_flo[2] ~ dnorm(0, 1/1000)

  #Colonisation
  beta_rho[1] ~ dnorm(0, 1/1000)
  beta_rho[2] ~ dnorm(0, 1/1000)
  
  #Survival
  beta_phi[1] ~ dnorm(0, 1/1000)
  beta_phi[2] ~ dnorm(0, 1/1000)
  
  #Colonization
  alpha ~ dgamma(0.001,0.001)   # was dunif(-10,10), Dispersal scale parameter
  #rho ~ dbeta(1,1)  # Population-level per capita effective dispersal rate

  ###Derived quantities block
  #First year
  psi[1] <- psi1
  n_occ[1]<- sum(z[1:n.sites,1])
  p_occ[1] <- n_occ[1]/n.sites
  #elev_max[1] <- max(elev[z[1:n.sites,1]>0])
  
  #All subsequent years
  for(t in 1:(t.max-1)) {
    n_occ[t+1] <- sum(z[1:n.sites,t+1])
    p_occ[t+1] <- n_occ[t+1]/n.sites
    #elev_max[t+1] <- max(elev[z[1:n.sites,t+1]])
      for(i in 1:n.sites){
        colyn[i,t+1] <- step(z[i,t+1]-z[i,t]-1)
        extyn[i,t+1] <- step(z[i,t]-z[i,t+1]-1)
        turyn[i,t+1] <- step(abs(z[i,t+1]-z[i,t])-1)
      }
    scol[t+1] <- sum(colyn[1:n.sites, t+1])
    sext[t+1] <- sum(extyn[1:n.sites,t+1])
    stur[t+1] <- sum(turyn[1:n.sites,t+1])
    }
  
  }
  ",
  fill = TRUE, file="occ1.txt")

sink()
#you could ask just for the mean, is in rjags! Look at means only, not posterior (example in WAIC folder in Day 3)
#Put in truncated dispersal kernel, estimate range limit


### 2) Set up a list that contains all the necessary data

Data_simple <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(jags_data$occ), t.max = ncol(jags_data$occ), y = jags_data$occ, N=jags_data$N, flo=jags_data$flo, elev=elev) 

# 3) Specify a function to generate inital values for the parameters
#inits_fn = function() list(alpha=0.1, rho=0.1, psi1 = 0.1, tau_dem= 1, tau_flo= 1, tau_pat = 1, beta_dem=runif(2,-3,3), beta_flo=runif(2,-3,3), beta_pat=runif(2,-3,3), beta_phi=runif(2,-3,3), p = 0.8) 
inits_fn = function() list(alpha=0.1, psi1 = 0.1, tau_N= 1, tau_flo= 1, beta_N=runif(2,-3,3), beta_flo=runif(2,-3,3), beta_phi=runif(2,-3,3), beta_rho=runif(2,-3,3))

jagsModel_nop = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 3, n.adapt= 0, inits = inits_fn)

# Specify parameters for which posterior samples are saved
para.names = c('n_occ', 'p_occ', 'alpha', 'scol', 'sext', 'stur',
               'beta_rho', 'beta_phi', 'beta_N', 'beta_flo') 

### 4) Continue the MCMC runs with sampling
Samples = coda.samples(jagsModel_nop, variable.names = para.names, n.iter = 10, na.rm=F)


par_vals = summary(Samples)$statistics
summary(Samples)
plot(Samples,ask=T)

save(Samples, file='modeloutput092020_solz.RData')

