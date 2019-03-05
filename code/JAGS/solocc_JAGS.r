library(rjags)

###########################################################
load('../data/soldat.RData')


# 2008, 2010, 2012
Occ <- t(stan_data$y)
Occ.dat <- matrix(NA, nrow = nrow(Occ), ncol = 5)
colnames(Occ.dat) <- (2008:2012)
Occ.dat[, c('2008','2010','2012')] <- Occ

D <- stan_data$D
diag(D) <- Inf

# define maximum distance over which neighbour patches are considered as potential 
# sources of colonization
max.dist <- 500

# for now, I kick out all patches that have no no neighbours within max.dist
# (for technical reasons)
rem.ind <- c(38,81,2697) 
Occ.dat <- Occ.dat[-rem.ind,]
D <- D[-rem.ind,-rem.ind]


n.sites <- nrow(Occ.dat)
t.max <- ncol(Occ.dat)

# Create a list of neighbours within max.dist for each patch
NB.list <- apply(D,1, function(x) which(x <= max.dist))
# Number of neighbours per patch
n.nb <- sapply(NB.list,length)
summary(n.nb)

# Here I create a matrix NB.mat with the patch numbers (indices) of all neighbours
# for each patch
# and a re-formated distance matrix D.nb with the distances to these neighbours 
NB.mat <- matrix(NA, nrow = n.sites, ncol = max(n.nb))
D.nb <- NB.mat
for(i in 1:n.sites){
  NB.mat[i,1:n.nb[i]] <- NB.list[[i]] 
  D.nb[i,1:n.nb[i]] <- D[i,NB.list[[i]]]
  }

###########################################################

# spatial model

# Save a description of the model in JAGS syntax to working directory
sink("SPOM.txt")
cat(
  "
  var S.nb[n.sites, t.max-1, max(n.nb)], C[n.sites, t.max-1], S[n.sites, t.max-1], occ[n.sites, t.max], p.occ[n.sites, t.max]
  
  model{

  # Data model
  for(t in 1:t.max){
    for(i in 1:n.sites){
      p.obs[i,t] <- occ[i,t] 
      y[i,t] ~ dbin(p.obs[i,t], 1)
    }
  }

  # Process model
  # 1st year 
  for(i in 1:n.sites){
    p.occ[i,1] <- psi0
    occ[i,1] ~ dbin(p.occ[i,1], 1)
    }
  
  
  for(t in 1:(t.max-1)){
  for(i in 1:n.sites){

  for(j in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
  S.nb[i,t,j] <- occ[NB.mat[i,j],t] * exp(-a * D.nb[i, j])
  }

  S[i,t] <- sum(S.nb[i,t,1:n.nb[i]])
  C[i,t] <- (S[i,t]*S[i,t]) / ( (S[i,t]*S[i,t]) + y2)
  
  p.occ[i,t+1] <- occ[i,t] * (1-E) + (1 - occ[i,t])*C[i,t] + 0.001
  occ[i,t+1] ~ dbin(p.occ[i,t+1], 1)
  }
  }


  # Prior distributions
  psi0 ~ dbeta(1,1)
  E ~ dbeta(1,1)
  a ~ dunif(0,5)
  y2 ~ dunif(0,100)
  }
  ",fill = TRUE)

sink()

###############

# 2) Set up a list that contains all the necessary data
Data = list(y = Occ.dat, D.nb = D.nb, NB.mat = NB.mat, n.nb = n.nb, n.sites = n.sites, t.max = t.max)

# 3) Specify a function to generate inital values for the parameters
# Initial values for occ are set to the value in Occ.dat, where data exixist,
# and to 1 otherwise
init.occ <- Occ.dat
init.occ[is.na(init.occ)] <- 1
inits.fn <- function() list(psi0 = 0.1, E = 0.1, y2 = 25, a = 1, occ = init.occ)

# Compile the model and run the MCMC for an adaptation (burn-in) phase
jagsModel <- jags.model(file= "SPOM.txt", data=Data, init = inits.fn, n.chains = 3, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names <- c("psi0","E","y2","a")

# Continue the MCMC runs with sampling
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

# Plot the mcmc chain and the posterior sample for p
plot(Samples)

# convergence check
coda::gelman.diag(Samples)

# Statistical summaries of the posterior sample
summary(Samples)

Samples.wd <- window(Samples, start = 1500 )
plot(Samples.wd)

#Plot predicted quantities over a range of years?
turnout.mat <- rbind(Samples[[1]], Samples[[2]], Samples[[3]])
Samp.out <- data.frame(turnout.mat)
p <- Samp.out$E

p.mean <- mean(p)

plot(p.mean ~ c(1:5), xlim = c(1,5))
