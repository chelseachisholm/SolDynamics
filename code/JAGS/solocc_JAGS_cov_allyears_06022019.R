library(rjags)
library(dplyr)
library(tidyr)

rm=list(ls())
###########################################################
load('./data/jags_data_allyears.RDS')

# 2008, 2009(NAs), 2010, 2011, 2012
Occ <- jags_data$C
Occ.dat <- matrix(NA, nrow = nrow(Occ), ncol = 5)
colnames(Occ.dat) <- c('2008','2009', '2010','2011','2012')
Occ.dat[, c('2008','2009', '2010','2011','2012')] <- Occ
Occ.dat[,2]<-NA

D <- as.matrix(jags_data$D)
diag(D) <- Inf

# define maximum distance over which neighbour patches are considered as potential 
# sources of colonization
max.dist <- 500

# for now, I kick out all patches that have no no neighbours within max.dist
# (for technical reasons)
rem.ind <- c(38,81,2697) 
Occ.dat <- Occ.dat[-rem.ind,]
Occ.dat[Occ.dat>1]<-1
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

#Create covariates
elev<- jags_data$elev
tann<- jags_data$tann
ddist<- jags_data$ddist
elev <- elev[-rem.ind,]
tann <- tann[-rem.ind,]
ddist <- ddist[-rem.ind,]
elev <- (elev-mean(elev)) / sd(elev)
tann <- (tann-mean(tann, na.rm=T)) / sd(tann, na.rm=T)
ddist <- (ddist-mean(ddist,na.rm=T)) / sd(ddist,na.rm=T)
logit <- function(x) {
  log(x/(1 - x))
}

###########################################################

# spatial model

# Save a description of the model in JAGS syntax to working directory
sink("SPOM.txt")
cat(
  "

  model{

  # Data model
  for(t in 1:t.max){
    for(i in 1:n.sites){
      muY[i,t] <- z[i,t] 
      y[i,t] ~ dbern(muY[i,t])
    }
  }

  # Process model
  # 1st year 
  for(i in 1:n.sites){
    z[i,1] ~ dbern(psi1)
  
  for(t in 1:(t.max-1)){
  
  #  Model of local survival (1-extinction) at site i
  logit(phi[i,t]) <- alpha.phi + beta.phi1*elev[i,t] 


  #  Pairwise 'source strength' calc and colonisation probability
      for(j in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
        S.nb[i,t,j] <- z[NB.mat[i,j],t] * exp(-a * D.nb[i, j])
      }

  S[i,t] <- sum(S.nb[i,t,1:n.nb[i]])
  C[i,t] <- (S[i,t]*S[i,t]) / ( (S[i,t]*S[i,t]) + y2)
  
  muZ[i,t+1] <- (z[i,t]) * (phi[i,t]) + (1 - z[i,t])*C[i,t] + 0.001
  z[i,t+1] ~ dbern(muZ[i,t+1])
    }
  }

 
 # Prior distributions
  psi1 ~ dbeta(1,1)
  alpha.phi ~ dnorm(0,0.3)
  beta.phi1 ~ dnorm(0, 0.3)
  a ~ dunif(0,5)
  y2 ~ dunif(0,100)
  # for (t in 1:(t.max-1)){
  # phi[t] ~ dunif(0, 1) 
  # C[t] ~ dunif(0, 1)
  # }
  
psi[1]<-psi1
n.occ[1]<-sum(z[1:n.sites,1])
for(t in 1:(t.max-1)) {
#  psi[t+1]<-psi[t]*phi[t] + (1-psi[t])*C[t] #issue because not site specific (make hyperparameter?)
  n.occ[t+1]<-sum(z[1:n.sites,t+1])
#  growthr[t+1]<-psi[t+1]/psi[t]
#  turnover[k]<-(1-psi[k])*C[k]/psi[k+1]
}
}
  ",fill = TRUE)

sink()

###############

# 2) Set up a list that contains all the necessary data
Data = list(y = Occ.dat, D.nb = D.nb, NB.mat = NB.mat, n.nb = n.nb, n.sites = n.sites, t.max = t.max, elev=elev)

# 3) Specify a function to generate inital values for the parameters
# Initial values for occ are set to the value in Occ.dat, where data exixist,
# and to 1 otherwise
init.occ <- Occ.dat
init.occ[is.na(init.occ)] <- 1
inits.fn <- function() list(psi1 = 0.1, alpha.phi=0.1, beta.phi1=0.1, phi=0.1, y2 = 25, a = 1, z = init.occ)

# Compile the model and run the MCMC for an adaptation (burn-in) phase
jagsModel <- jags.model(file= "SPOM.txt", data=Data, n.chains = 3, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names <- c("alpha.phi","beta.phi1","y2","a", "n.occ") #could monitor phi or z but will save each sample from each site (so long!)
  #extract averages across bins of elevation? so for every bin, do this, something like that, so we have estimates for each bin? Or is there a way to just simulate directly
  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

# Continue the MCMC runs with sampling
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

save.image('jags_ouput_allyears.RData')

#### TO INCLUDE IN FUTURE ####

#A derived quantities block (but per elevation band??) Ex:
#Derived quantities
psi[0] <- psi0 
n.occ[1]<-sum(occ[1:nsite,1]) 
for (t in 2:t.max){
  psi[t] <- (occ[t]) * (phi[t]) + (1 - occ[t])*C[t] + 0.001 
  n.occ[t] <- sum(occ[1:nsite,t])
  growthr[t] <- phi[t]/phi[t−1]
  turnover[t−1] <- (1 − phi[t-1]) * gamma[t−1]/phi[t]
} }
