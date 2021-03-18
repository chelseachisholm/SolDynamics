#### Full model 30.03.2020
library(rjags)
# Load data
load('./jags_data.RData')

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

#Remove 4 cells with 0 values in N and flo (check these?)
jags_data$N[jags_data$N==0]<-NA
jags_data$flo[jags_data$flo==0]<-NA

#transform N and flo
# jags_data$N <- log10(jags_data$N)
# jags_data$flo <- log10(jags_data$flo)

#Create NA columns to simulate data to 2020 (+8)
# jags_data$occ <- as.matrix(cbind(jags_data$occ, matrix(NA, 2741, 8)))
# jags_data$z <- as.matrix(cbind(jags_data$z, matrix(NA, 2741, 8)))
# jags_data$N <- as.matrix(cbind(jags_data$N, matrix(NA, 2741, 8)))
# jags_data$flo <- as.matrix(cbind(jags_data$flo, matrix(NA, 2741, 8)))
# jags_data$nyear <- 13

#simulations
#reduce N across all
#jags_data$N[,5] <- jags_data$N[,5]/2
#reduce pop number
#bo <- sample(jags_data$occ(which(jags_data$occ[,5]==1)), 112) #224 total, so 112 is half
#jags_data$occ[,5]

#load.module('glm')

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
  
  # Modelling missing data for demographic data (data scaled!--> zeros so could use lognormal to constrain to positive values)
  N[i,t] ~ dnorm(mu_N[i,t], tau_N) 
  flo[i,t] ~ dnorm(mu_flo[i,t], tau_flo) 

  # Models of elevation effects on dem data
  log(mu_N[i,t]) <- beta_N[1] + beta_N[2]*elev[i] 
  log(mu_flo[i,t]) <- beta_flo[1] + beta_flo[2]*elev[i] 
  
  # Creating neighboring patch weights by total seed output per patch
  logit(totseed[i,t]) <- N[i,t]*flo[i,t]
  
  # Neighbourhood connectivity to focal patch i (weighted by total seed production of neighbourhood)
    for(n in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
      gammaDistPairs[i,n,t] <- totseed[NB.mat[i,n],t] * exp(-alpha*D.nb[i,n]) * z[NB.mat[i,n],t]
    }
  
  # Source strength calc
  S[i,t] = prod(1-gammaDistPairs[i,1:n.nb[i],t]) #site-level 'connectivity'

  # Model of colonization probability 
  gamma[i,t] = 1 - exp(-rho[i,t]*S[i,t])

  logit(rho[i,t]) <- beta_rho[1] + beta_rho[2]*elev[i] + trand[t] 

  # Model of local survival probability (1-extinction) + random year effect 
  logit(phi[i,t]) <- beta_phi[1] + beta_phi[2]*elev[i] + trand[t] 
  
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
  taut ~ dgamma(0.001,0.001)
  
  #For missing data in flowering and dem data
  sigma_N ~ dunif(0.00001, 10000) #N is highly variable
  tau_N <- 1/(sigma_N*sigma_N)
  sigma_flo ~ dunif(0.00001, 100)
  tau_flo <- 1/(sigma_flo*sigma_flo)
  
  #Demography
  beta_N[1] ~ dnorm(0, 0.137)
  beta_N[2] ~ dnorm(0, 0.137)
  beta_flo[1] ~ dnorm(0, 0.137)
  beta_flo[2] ~ dnorm(0, 0.137)

  #Colonisation
  beta_rho[1] ~ dnorm(0, 0.001)
  beta_rho[2] ~ dnorm(0, 0.001)
  
  #Survival
  beta_phi[1] ~ dnorm(0, 0.001)
  beta_phi[2] ~ dnorm(0, 0.001)
  
  #Colonization
  log.alpha ~ dnorm(0, 0.137)   #dispersal scale parameter
  alpha <- exp(log.alpha)

  ###Derived quantities block
  #First year
  psi[1] <- psi1
  n_occ[1]<- sum(z[1:n.sites,1])
  p_occ[1] <- n_occ[1]/n.sites
  
  #All subsequent years
  for(t in 1:(t.max-1)) {
    n_occ[t+1] <- sum(z[1:n.sites,t+1])
    p_occ[t+1] <- n_occ[t+1]/n.sites
    # #elev_max[t+1] <- max(elev[z[1:n.sites,t+1]])
    #   for(i in 1:n.sites){
    #     colyn[i,t+1] <- step(z[i,t+1]-z[i,t]-1)
    #     extyn[i,t+1] <- step(z[i,t]-z[i,t+1]-1)
    #     turyn[i,t+1] <- step(abs(z[i,t+1]-z[i,t])-1)
    #   }
    # scol[t+1] <- sum(colyn[1:n.sites, t+1])
    # sext[t+1] <- sum(extyn[1:n.sites,t+1])
    # stur[t+1] <- sum(turyn[1:n.sites,t+1])
    }
  
  }
  ",
  fill = TRUE, file="occ1.txt")

sink()
#you could ask just for the mean, is in rjags! Look at means only, not posterior (example in WAIC folder in Day 3)
#Put in truncated dispersal kernel, estimate range limit


### 2) Set up a list that contains all the necessary data

Data_simple <- list(n.nb = n.nb, NB.mat = NB.mat, D.nb = D.nb, n.sites = nrow(jags_data$occ), t.max = ncol(jags_data$occ), y = jags_data$occ, N=jags_data$N, flo=jags_data$flo, elev=elev) 

# Specify a function to generate inital values for the parameters
inits_fn = function() list(z=matrix(data=rep(1, 2741*5), nrow=2741, ncol=5), log.alpha=0.001, psi1 = 1, taut = 0.001, sigma_N= runif(1, 0, 5), sigma_flo= runif(1, 0, 5), beta_N=rnorm(2,0,0.001), beta_flo=runif(2,0,0.001), beta_phi=runif(2,0,0.137), beta_rho=runif(2,0,0.137))

# Specify parameters for which posterior samples are saved
para.names = c('n_occ', 'p_occ', 'log.alpha', 'sigma_N', 'sigma_flo', 'taut',
               'beta_rho', 'beta_phi', 'beta_N', 'beta_flo') 

library('runjags')
parsamples <- run.jags('occ1.txt', data=Data_simple, inits = inits_fn, n.chains = 3, monitor=para.names,
                       adapt = 1000, burnin = 4000, sample = 5000, method='rjparallel', jags.refresh=10)
summary(parsamples)
newparsamples <- extend.jags(parsamples, sample=5000)
summary(newparsamples)
plot(newparsamples)
parsamples$mcmc
plot(parsamples$mcmc)
parsamples$summaries
#extract(parsamples, what = "dic")
#extract(parsamples, what = "ped")
gelman.diag(parsamples$mcmc)


# #using rjags
# jagsModel_nop_v2 = jags.model(file= "occ1.txt", data=Data_simple, n.chains = 2, inits = inits_fn)
# update(jagsModel_nop_v2, 500)
# Samples = coda.samples(jagsModel_nop_v2, variable.names = para.names, n.iter = 1000)
# par_vals = summary(Samples)$statistics
# summary(Samples)
# plot(Samples,ask=T)
# plot(Samples)


save(newparsamples, file='modeloutput012021.RData')

