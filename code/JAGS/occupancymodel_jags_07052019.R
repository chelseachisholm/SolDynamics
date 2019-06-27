#### RUN MODEL ####

Run_JAGSmodel <- function(x) {

# 2008, 2009(NAs), 2010, 2011, 2012
Occ <- x$C
Occ.dat <- matrix(NA, nrow = nrow(Occ), ncol = 5)
colnames(Occ.dat) <- c('2008','2009', '2010','2011','2012')
Occ.dat[, c('2008','2009', '2010','2011','2012')] <- Occ
Occ.dat[,2]<-NA

D <- as.matrix(x$D)
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
elev<- x$elev
tann<- x$tann
ddist<- x$ddist
elev <- elev[-rem.ind,]
#tann <- tann[-rem.ind,]
ddist <- ddist[-rem.ind,]
elev <- (elev-mean(elev)) / sd(elev)
#tann <- (tann-mean(tann, na.rm=T)) / sd(tann, na.rm=T)
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

  # Observation model
  for(t in 1:t.max){
    for(i in 1:n.sites){
      muY[i,t] <- z[i,t] 
      y[i,t] ~ dbern(muY[i,t])
    }
  }

  # State (process) model
  # 1st year 
  for(i in 1:n.sites){
    z[i,1] ~ dbern(psi1)
  
  for(t in 1:(t.max-1)){
  
  #  Model of local survival (1-extinction) at site i
  logit(phi[i,t]) <- beta.phi[1] + beta.phi[2] * elev[i,t] 


  #  Pairwise 'source strength' calc and colonisation probability
      for(j in 1:n.nb[i]){ # loop of the nb[i] neighbours of patch i
        S.nb[i,t,j] <- z[NB.mat[i,j],t] * exp(-a * D.nb[i, j])
      }

  S[i,t] <- sum(S.nb[i,t,1:n.nb[i]])
  C[i,t] <- (S[i,t]*S[i,t]) / ( (S[i,t]*S[i,t]) + y2)
  
  #defining occupancy probability
  muZ[i,t+1] <- (z[i,t]) * (phi[i,t]) + (1 - z[i,t])*C[i,t] + 0.001
  #new occupancy probability  
  z[i,t+1] ~ dbern(muZ[i,t+1])
    }
  }

 
 # Prior distributions
  psi1 ~ dbeta(1,1)
  beta.phi[1] ~ dnorm(0,0.3)
  beta.phi[2] ~ dnorm(0, 0.3)
  a ~ dunif(0,5)
  y2 ~ dunif(0,100)
  # for (t in 1:(t.max-1)){
  # phi[t] ~ dunif(0, 1) 
  # C[t] ~ dunif(0, 1)
  # }
  
#Derived quantities block
psi[1]<-psi1
n.occ[1]<-sum(z[1:n.sites,1])
for(t in 1:(t.max-1)) {
    n.occ[t+1]<-sum(z[1:n.sites,t+1])
  #   for(i in 1:n.sites) {
  # psi[i,t+1]<-z[i,t]*phi[i,t] + (1-z[i,t])*C[i,t] #issue because not site specific (make hyperparameter?)
  # growthr[i,t+1]<-psi[t+1]/psi[t]
  # turnover[i,t]<-(1-psi[i,t])*C[i,t]/psi[i,t+1]
  #  }#i
  }#t
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
inits.fn <- function() list(psi1 = 0.1, beta.phi=runif(2,-3,3), phi=0.1, y2 = 25, a = 1, z = init.occ)

# Compile the model and run the MCMC for an adaptation (burn-in) phase
jagsModel <- jags.model(file= "SPOM.txt", data=Data, n.chains = 3, n.adapt= 1000)
# Specify parameters for which posterior samples are saved
para.names <- c("beta.phi","y2","a", "n.occ") #could monitor phi or z but will save each sample from each site (so long!)
  #extract averages across bins of elevation? so for every bin, do this, something like that, so we have estimates for each bin? Or is there a way to just simulate directly
  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...

# Continue the MCMC runs with sampling
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)

return(Samples)
}

