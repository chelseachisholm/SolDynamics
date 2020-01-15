#### RUN MODEL ####

Run_JAGSmodel <- function(x) {
  
  
  source("./helperfunctions.R")
  
  # 2008, 2009(NAs), 2010, 2011, 2012
  x<-jags_data
  occ <- as.matrix(x$occ)
  nNeighbours<- x$nsite
  n.sites <- x$nsite
  t.max <- ncol(occ)
  dist.all <- as.matrix(x$D)
  dem<-x$dem
  dem[,c(1,2,4)]<- dem[,3]
  flo<-x$flo
  
 dem[is.na(dem)] <- median(dem, na.rm=T)
 flo[is.na(flo)] <- median(flo, na.rm=T)
 #flo[is.na(flo)] <- runif(length(is.na(flo)), 0,10)
  
  #Create covariates
  elev<- zscale(x$elev)
  tann<- zscale(x$tann)
  
  #Create list of neighbours within max distance (using 500 for now)
  # Create a list of neighbours within max.dist for each patch
  # max.dist <- 30000
  # NB.list <- apply(x$D,1, function(x) which(x <= max.dist))
  # # Number of neighbours per patch
  # n.nb <- sapply(NB.list,length)
  # summary(n.nb)
  # 
  # # Here I create a matrix NB.mat with the patch numbers (indices) of all neighbours
  # # for each patch
  # # and a re-formated distance matrix D.nb with the distances to these neighbours 
  # NB.mat <- matrix(NA, nrow = n.sites, ncol = max(n.nb))
  # D.nb <- NB.mat
  # for(i in 1:n.sites){
  #   NB.mat[i,1:n.nb[i]] <- NB.list[[i]] 
  #   D.nb[i,1:n.nb[i]] <- D[i,NB.list[[i]]]
  # }
  # 
  logit <- function(x) {
    log(x/(1 - x))
  }
  
  ###########################################################
  
  # spatial model (old)
  
  # Save a description of the model in JAGS syntax to working directory
  sink("occ1.txt")
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
          logit(phi[i,t]) <- beta_phi[1] + beta_phi[2] * elev[i,t] #+ beta_phi[3] * tann[i,t]
    
    #Pairwise colonization probability based on dispersal kernel 
		for(n in 1:nNeighbours){ 
    #Exponential kernel, probability of dispersal to a site based on distance
      gammaDistPairs[i,n,t] = gamma0 * exp(-gamma0*dist.all[i,n]) * z[n,t] 
    }#add in repro ability for other patches

    # Colonization prob
    gamma[i,t] = 1 - prod(1-gammaDistPairs[i,,t])

    ##Flowering:
    logit(kappa[i,t]) = beta_f[1] + beta_f[2] * flo[i,t] 

    ##Relative abundance
    logit(lambda[i,t]) = beta_a[1] + beta_a[2] * dem[i,t]
    
    Ez[i,t+1] = kappa[i,t]*lambda[i,t]*gamma[i,t]*(1-z[i,t]) + phi[i,t]*z[i,t]+0.00001
		#This won't normalize without the small additive constant. Why? 
    

    #new occupancy probability  
    z[i,t+1] ~ dbern(Ez[i,t+1])
    }
    }
    
    
    # Prior distributions
	  psi1 ~ dbeta(1,1)

    #priors for missing data in predictors
    #dem[i,t] ~ dnorm(0, 0.03)  
    #flo[i,t] ~ dnorm(0, 0.03) 


    #For predictors
    ###Survival
    beta_phi[1] ~ dnorm(0,0.3)
    beta_phi[2] ~ dnorm(0, 0.3)
    #beta_phi[3] ~ dnorm(0, 0.3)
    
    ###Flowering
    beta_f[1] ~ dnorm(0,0.3)
    beta_f[2] ~ dnorm(0, 0.3)
    
    ###Relative abundance
    beta_a[1] ~ dnorm(0,0.3)
    beta_a[2] ~ dnorm(0, 0.3)

    ###Dispersal
    #gamma0 ~ dbeta(1,1)
    gamma0 ~ dnorm(0,0.1)  
    
    ###Derived quantities block
    psi[1]=psi1
    n_occ[1]=sum(z[1:n.sites,1])
    
    for(t in 1:(t.max-1)) {
    n_occ[t+1]=sum(z[1:n.sites,t+1])
#also add range limit
    
	}
}
",fill = TRUE, file="occ1.txt")
  
  sink()
  
  #Put in truncated dispersal kernel, estimate range limit
  
  
  
  ###############
  
  ### 2) Set up a list that contains all the necessary data
  Data = list(y = occ, dist.all = dist.all,  n.sites = n.sites, elev = elev, #tann = tann,  
              nNeighbours=n.sites, t.max = t.max, flo = flo, dem = dem)
  
  # 3) Specify a function to generate inital values for the parameters
  init.pa = occ
  inits_fn = function() list(psi1 = 0.1, beta_phi=runif(2,-3,3), beta_f=runif(2,-3, 3), 
                             beta_a=runif(2,-3, 3), gamma0 = 0.1, z = init.pa)
  #inits_fn = function() list(psi1 = 0.1, beta_phi=runif(1,-3,3), beta_f=runif(3,-3,3), z = init.pa)
  load.module('glm')
  jagsModel = jags.model(file= "occ1.txt", data=Data, n.chains = 2, n.adapt= 1000)
  # Specify parameters for which posterior samples are saved
  para.names = c("beta_phi","beta_f","beta_a","psi1","gamma0", "n_occ")  #all the data for one parameter of interest, like colonization probability, using some of the other parameter estimes. Hmmm...
  
  ### 4) Continue the MCMC runs with sampling
  Samples = coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)
  
  ### 5) Check the outpt values, and look for convergence:
  par_vals = summary(Samples)$statistics
  summary(Samples)
  plot(Samples,ask=T)
  
  ### 6) Compare to values in the underlying population model:
  p_surv=antilogit(par_vals[5,1]+par_vals[6,1]*elev)
  p_birth=antilogit(par_vals[3,1]+par_vals[4,1]*flo)
  p_comp = antilogit(par_vals[1,1]+par_vals[2,1]*dem)
  
  plot(elev, p_surv)
  plot(p_birth, elev)
  plot(p_comp, elev)
  
  
  
  #An approximation of competition from model: 
  #comp2 =1/(1+nr_act[,,2]*alphas[1,2])
  #Total expected competitive effect from occupancy model (prob_flowers*prob_comp): 
  p_comp_tot = p_birth*p_comp
  
  #Plot
  plot(p_comp_tot[,1])
  points(p_comp[,1],col="red")
  
  #Forecast (modify this code!)
  jags.data = list("y"=c(occ,NA,NA,NA),"t.max"=(t.max+3))
  jags.params=c("sd.q","sd.r","predY","mu")
  model.loc=("ss_model.txt")
  mod_ss_forecast = jags(jags.data, parameters.to.save=jags.params,
                         model.file=model.loc, n.chains = 3, n.burnin=5000, n.thin=1,
                         n.iter=10000, DIC=TRUE)
  
  
  return(Samples)
}

