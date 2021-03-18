
#### Lognormal distribution without zeros ####
dat <- jags_data$N[,3]
sink("occ1.txt") 

cat(
  "
  model{
  for(i in 1:N){ # loop through all data points
  y[i] ~ dnorm(mu[i], tau)
  log(mu[i]) <-beta.l[1] + beta.l[2]*elev[i] 
  
  }
  # priors:
  beta.l[1] ~ dnorm(0, 0.001) # overall model intercept
  beta.l[2] ~ dnorm(0, 0.001) #slope
  tau <- pow(sigma, -2)                                          # precision
  sigma ~ dunif(0, 10000)                                          # standard deviation
  }",
fill = TRUE, file="occ1.txt")

sink()

parameters <- c("beta.l", "sigma") # which parameters are we interested in getting reported?
inits <- function() list(beta.l=runif(2, 0, 2), sigma=2)
dat.l <- list(N=2741, y=dat, elev=elev)
LRjags <- run.jags('occ1.txt', data=dat.l, inits=inits, n.chains = 3, monitor=parameters,
                   burnin = 1000, sample = 4000, method='rjparallel')
summary(LRjags)
plot(LRjags)

#### Log normal with new priors ####
dat <- log(jags_data$flo[,3])
sink("occ1.txt") 

cat(
  "
  model{
  for(i in 1:N){ # loop through all data points
  y[i] ~ dnorm(mu, tau)
  
  }
  # priors:
  tau <- pow(sigma, -2)                                          # precision
  sigma ~ dexp(0.5) #0.5 for longer tail  
  mu ~ dnorm(0, 0.1)
  }",
fill = TRUE, file="occ1.txt")

sink()

parameters <- c("mu", "sigma") # which parameters are we interested in getting reported?
inits <- function() list(mu=rnorm(1, 0, 0.1), sigma=2)
dat.l <- list(N=2741, y=dat)
LRjags <- run.jags('occ1.txt', data=dat.l, inits=inits, n.chains = 3, monitor=parameters,
                   burnin = 1000, sample = 4000, method='rjparallel')
summary(LRjags)
plot(LRjags)

#Prepping data, checking frequentist approaches
dat <- jags_data$N[,3]
dat[is.na(dat)]<-0
m1<-glm(dat~elev,family="poisson")
summary(m1)
plot(m1)

#Checking distributions of data
N <- jags_data$N[,3]
hist(N)
flo <- jags_data$flo[,3]
hist(flo)

#### Zero inflated poisson distribution ####
sink("occ1.txt") 

cat(
  "
  model{
  for(i in 1:N){ # loop through all data points
    y[i] ~ dpois(mu[i])
    mu[i] <- lambda[i]*z[i] #+ 0.00001 ## hack required for Rjags -- otherwise 'incompatible'-error
    z[i] ~ dbern(psi)
    log(lambda[i]) <-beta.l[1] + beta.l[2]*elev[i] 

  }
  # priors:
  beta.l[1] ~ dnorm(0, 0.01) # overall model intercept
  beta.l[2] ~ dnorm(0, 0.01) #slope
  psi ~ dunif(0, 1) # proportion of non-zeros
}",
fill = TRUE, file="occ1.txt")

sink()

parameters <- c("beta.l", "psi") # which parameters are we interested in getting reported?
inits <- function() list(beta.l=runif(2, 0, 2), psi = runif(1, 0, 1))
dat.l <- list(N=2741, y=dat, elev=elev, z = jags_data$occ[,3])
ZIPRjags <- run.jags('occ1.txt', data=dat.l, inits=inits, n.chains = 3, monitor=parameters,
                    burnin = 1000, sample = 4000, method='rjparallel')
summary(ZIPRjags)

plot(ZIPRjags)

#### Poisson distribution ####
sink("occ1.txt") 

cat(
  "
  model{
  for(i in 1:N){ # loop through all data points
  y[i] ~ dpois(mu[i])
  mu[i] <- lambda[i] # + 0.00001 ## hack required for Rjags -- otherwise 'incompatible'-error

  log(lambda[i]) <-beta.l[1] + beta.l[2]*elev[i] 
  
  }
  # priors:
  beta.l[1] ~ dnorm(0, 0.01) # overall model intercept
  beta.l[2] ~ dnorm(0, 0.01) #slope
  }",
fill = TRUE, file="occ1.txt")

sink()

#with zeros
parameters <- c("beta.l") # which parameters are we interested in getting reported?
inits <- function() list(beta.l=runif(2, 0, 2))
dat.l <- list(N=2741, y=dat, elev=elev)
PRjags <- run.jags('occ1.txt', data=dat.l, inits=inits, n.chains = 3, monitor=parameters,
                     burnin = 1000, sample = 4000, method='rjparallel')
summary(PRjags)

plot(PRjags)

#### Poisson distribution without zeros in data ####
#removing zeros
parameters <- c("beta.l") # which parameters are we interested in getting reported?
inits <- function() list(beta.l=runif(2, 0, 2))
y.n  <- dat[jags_data$occ[,3]>0]
elev.n <-elev[jags_data$occ[,3]>0]
length(y)

dat.l <- list(N=123, y=y.n, elev=elev.n)
PRjags <- run.jags('occ1.txt', data=dat.l, inits=inits, n.chains = 3, monitor=parameters,
                   burnin = 1000, sample = 4000, method='rjparallel')
summary(PRjags)

plot(PRjags)
m1<-glm(y.n~elev.n,family="poisson")
summary(m1)

#### Lognormal distribution ####

sink("occ1.txt") 

cat(
  "
  model{
  for(i in 1:N){ # loop through all data points
  y[i] ~ dnorm(mu[i], tau)
  log(mu[i]) <-beta.l[1] + beta.l[2]*elev[i] 
  
  }
  # priors:
  beta.l[1] ~ dnorm(0, 0.001) # overall model intercept
  beta.l[2] ~ dnorm(0, 0.001) #slope
  tau <- pow(sigma, -2)                                          # precision
  sigma ~ dunif(0, 10000)                                          # standard deviation
  }",
fill = TRUE, file="occ1.txt")

sink()

parameters <- c("beta.l", "sigma") # which parameters are we interested in getting reported?
inits <- function() list(beta.l=runif(2, 0, 2), sigma=2)
dat.l <- list(N=2741, y=dat, elev=elev)
LRjags <- run.jags('occ1.txt', data=dat.l, inits=inits, n.chains = 3, monitor=parameters,
                     burnin = 1000, sample = 4000, method='rjparallel')
summary(LRjags)
plot(LRjags)

#### Zero-inflated lognormal dist ####
sink("occ1.txt") 

cat(
  "
  model{
  for(i in 1:N){ # loop through all data points
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- lambda[i]*z[i] #+ 0.00001 ## hack required for Rjags -- otherwise 'incompatible'-error
  z[i] ~ dbern(psi)
  log(lambda[i]) <-beta.l[1] + beta.l[2]*elev[i] 
  
  }
  # priors:
  beta.l[1] ~ dnorm(0, 0.001) # overall model intercept
  beta.l[2] ~ dnorm(0, 0.001) #slope
  tau <- pow(sigma, -2)                                          # precision
  sigma ~ dunif(0, 10000)     
  psi ~ dunif(0, 1) # proportion of non-zeros
  }",
fill = TRUE, file="occ1.txt")

sink()

parameters <- c("beta.l", "psi", "sigma") # which parameters are we interested in getting reported?
inits <- function() list(beta.l=runif(2, 0, 2), psi = runif(1, 0, 1))
dat.l <- list(N=2741, y=dat, elev=elev, z = jags_data$occ[,3])
ZIPRjags <- run.jags('occ1.txt', data=dat.l, inits=inits, n.chains = 3, monitor=parameters,
                     burnin = 2000, sample = 5000, method='rjparallel')
summary(ZIPRjags)

plot(ZIPRjags)


#### Two-stage lognormal dist ####
sink("occ1.txt") 

cat(
  "
  model{
  for(i in 1:N){ # loop through all data points
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- lambda[i]*z[i] 
  z[i] ~ dbern(psi)
  log(lambda[i]) <-beta.l[1] + beta.l[2]*elev[i] 
  
  }
  # priors:
  beta.l[1] ~ dnorm(0, 0.001) # overall model intercept
  beta.l[2] ~ dnorm(0, 0.001) #slope
  tau <- pow(sigma, -2)                                          # precision
  sigma ~ dunif(0, 10000)     
  psi ~ dunif(0, 1) # proportion of non-zeros
  }",
fill = TRUE, file="occ1.txt")

sink()

parameters <- c("beta.l", "psi", "sigma") # which parameters are we interested in getting reported?
inits <- function() list(beta.l=runif(2, 0, 2), psi = runif(1, 0, 1))
dat.l <- list(N=2741, y=dat, elev=elev, z = jags_data$occ[,3])
ZIPRjags <- run.jags('occ1.txt', data=dat.l, inits=inits, n.chains = 3, monitor=parameters,
                     burnin = 2000, thin=5, sample = 10000, method='rjparallel')
summary(ZIPRjags)

plot(ZIPRjags)

#### Two-stage gamma dist ####
sink("occ1.txt") 

cat(
  "
  model{
  for(i in 1:N){ # loop through all data points
  #y[z==0] <- NA
  z[i] ~ dbern(psi)
  if (z[i] > 0) {
  y[i] ~ dnorm(mu[i], tau)
  log(mu[i]) <-beta.l[1] + beta.l[2]*elev[i]
  }
  
  }
  # priors:
  beta.l[1] ~ dnorm(0, 0.001) # overall model intercept
  beta.l[2] ~ dnorm(0, 0.001) #slope
  tau <- pow(sigma, -2)                                          # precision
  sigma ~ dunif(0, 10000)     
  psi ~ dunif(0, 1) # proportion of non-zeros
  }",
fill = TRUE, file="occ1.txt")

sink()

parameters <- c("beta.l", "psi", "sigma") # which parameters are we interested in getting reported?
inits <- function() list(beta.l=runif(2, 0, 2), psi = runif(1, 0, 1))
dat.l <- list(N=2741, y=dat, elev=elev, z = jags_data$occ[,3])
ZIPRjags <- run.jags('occ1.txt', data=dat.l, inits=inits, n.chains = 3, monitor=parameters,
                     burnin = 2000, thin=5, sample = 10000, method='rjparallel')
summary(ZIPRjags)

plot(ZIPRjags)

