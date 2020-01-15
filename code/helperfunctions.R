#=============================================================================
# This code contains basic functions to calculate species intrinsic ranges
# and determine components of the invasion growth rates for species.
#
# The model is assumed to be the spatially explicit Leslie-Gower model, where
# the intrinsic rate of reproduction is assumed to be the spatially-varying
# parameter. 
#
# Dispersal kernels are allowed to vary in space and be location-specific
#
# The functions are broken into 3 sections: 
# 1. Miscellanerous functions: Mostly helper functions for later routines.
#		Also some functions to set up theoretical intrinsic ranges for 
#		simulations. 
# 2. Invader and resident equilibrium densities: Functions to calculate 
#		components of the Low-density growth rates (LGR) in various ways
#		(including both analytical and numerical approaches).
#		Two subsections: 1) Multispecies (full community) case and 2) Pairwise
#		case. 
# 3. Environmental distance and site-LGR impact: Functions to determine 1)
#		1) the environmental distance between a site's current and future state
#		(environmental distance is the Euclidean distance in the multivariate
#		abiotic space), and its environmental/physical distance from the 
#		nearest analogue and, 
#		2) the incluence that each site has on the various LGR metrics of (2)
#
#=============================================================================


#=============================================================================
# FUNCTION DEFINITIONS
#=============================================================================

#=============================================================================
#Miscellaneous functions
#=============================================================================

#R equivalent for matlab's meshgrid
meshgrid=function(a,b) {
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
} 

#Logit function
logit <- function(x) {
  log(x/(1 - x))
}

#Inverse logit: 
antilogit <- function(x) {
  exp(x)/(1 + exp(x))
}

#Scale data by mean and SD: 
zscale = function(x) {
  (x - mean(x,na.rm=T))/ sd(x,na.rm=T)
}

#These two functions can be used to find upper and lower edges
#of the range distributions. Needs the peak position, variance
#of the distribution, and a tolerance specifying number
#of standard deviations to include. 
get.upper = function (np, burns, ngens, pks, Ds, D.tol= 3) {
  
  ngenst=burns+ngens
  #Test if the peak and spatial width parameters are a single value
  #(i.e. stationary or moving)
  if(length(pks)>1){
    pk1 = pks[1]
    pk.end = pks[2]
    pk.by = pks[3]
    peak.stc.burn=c(matrix(pk1, ngens,1))
    peak.stc.gens=seq(pk1, pk.end, by=pk.by)[1:burns]
    peak.stc = c(peak.stc.burn, peak.stc.gens)} else {
      pk1 = pks
      peak.stc=c(matrix(pk1, ngenst,1))
    }
  
  
  if(length(Ds)>1){
    Ds1 = Ds[1]
    Ds.end = Ds[2]
    Ds.by = Ds[3]
    Ds.stc.burn=c(matrix(Ds1, ngens,1))
    Ds.stc.gens=seq(Ds1, Ds.end, by=Ds.by)[1:ngens]
    Ds.stc= c(Ds.stc.burn, Ds.stc.gens)
  } else { 
    
    Ds1 = Ds
    Ds.stc=c(matrix(Ds1, ngens,1))
  } 
  
  sds.up=round(np/2)+peak.stc+sqrt(Ds.stc)*D.tol
  return(sds.up)
  
}

get.lower = function (np, burns, ngens, pks, Ds, D.tol= 3) {
  
  ngenst=burns+ngens
  #Test if the peak and spatial width parameters are a single value
  #(i.e. stationary or moving)
  if(length(pks)>1){
    pk1 = pks[1]
    pk.end = pks[2]
    pk.by = pks[3]
    peak.stc.burn=c(matrix(pk1, ngens,1))
    peak.stc.gens=seq(pk1, pk.end, by=pk.by)[1:burns]
    peak.stc = c(peak.stc.burn, peak.stc.gens)} else {
      pk1 = pks
      peak.stc=c(matrix(pk1, ngenst,1))
    }
  
  
  if(length(Ds)>1){
    Ds1 = Ds[1]
    Ds.end = Ds[2]
    Ds.by = Ds[3]
    Ds.stc.burn=c(matrix(Ds1, ngens,1))
    Ds.stc.gens=seq(Ds1, Ds.end, by=Ds.by)[1:ngens]
    Ds.stc= c(Ds.stc.burn, Ds.stc.gens)
    
  } else { 
    
    Ds1 = Ds
    Ds.stc=c(matrix(Ds1, ngens,1))
  } 
  
  sds.low=round(np/2)+peak.stc-sqrt(Ds.stc)*D.tol
  return(sds.low)
  
}


#=============================================================================
#Ranges and range shifts 
#=============================================================================
#Make a species' intrinsic range using a Gaussian distribution. 
#Variables correspond to: 
#	Fr	Max height of distribution (max reproduction)
#	pks	Vector with: 
#		 1. Initial spatial mean (position)
#		 2. pk.ends: Ending spatial mean
#		 3. by: Rate at which mean moves 
#		If not a vector, it is assumed stationary (no change in position) 
#
#	Ds	Vector with:
#		 1. Spatial variance (width)
#		 2. Ds.end: ending spatial variance
#		 3. Rate at which variance changes
#		If not a vector, it is assumed stationary (no change in position) 
#
#	stc	spatial array with space and time coordinates
#	burns	time period where distribution is stationary
#	ngens	time period of change for distribution
#=============================================================================

make.range=function(Fr, pks, Ds, stc, ngens, burns, fast.by = FALSE, dnwidth=NULL) { 
  
  
  np=dim(stc)[2]
  nt=dim(stc)[1]
  Fr.tmp = matrix(0, nt, np)
  
  #Test if the peak and spatial width parameters are a single value
  if(length(pks)>1){pk1 = pks[1]} else { 
    pk1 = pks} 
  
  if(length(Ds)>1){Ds1 = Ds[1]} else { 
    Ds1 = Ds} 
  
  #Burn phase
  X1=stc[1:burns,(1:(np)),2]
  G1D=exp(-((X1-pk1)^2)/(2*Ds1))
  G1D[G1D<1e-7]=0
  G1D=G1D/(max(G1D))*Fr
  Fr.tmp[1:burns,]=G1D
  
  #Test if it is the fast or normal scenario
  # Note: fast scenario inherently assumes changing peak position
  # BUT NO change in width. 
  
  if(fast.by==FALSE) {
    #Moving phase
    #Test if mean and variance change and create appropriate variables if
    #they do
    if(length(pks)>1){ pk.end = pks[2]
    pk.by = pks[3]
    peak.stc=seq(pk1, pk.end, by=pk.by)[1:ngens]} else {
      peak.stc=c(matrix(pk1, ngens,1))}
    
    
    if(length(Ds)>1){  Ds.end = Ds[2]
    Ds.by = Ds[3]
    Ds.stc=seq(Ds1, Ds.end, by=Ds.by)[1:ngens]} else {
      Ds.stc=c(matrix(Ds1, ngens,1))}
    
    
    X2=stc[(1+burns):nt,(1:(np)),2]
    G1D=exp(-((X2-peak.stc)^2)/(2*Ds.stc))
    G1D[G1D<1e-7]=0
    G1D=G1D/(max(G1D))*Fr
    Fr.tmp[(1+burns):nt,]=G1D
  } else { 
    
    pk.end = pks[2]
    pk.by = pks[3]
    peak.stc1=seq(pk1, pk.end, by=pk.by)[1:dnwidth]
    peak.stc2=matrix(pk.end, (ngens-dnwidth),1)
    peak.stc=c(peak.stc1,peak.stc2)
    
    Ds.stc=c(matrix(Ds1, ngens,1))
    
    X2=stc[(1+burns):nt,(1:(np)),2]
    G1D=exp(-((X2-peak.stc)^2)/(2*Ds.stc))
    G1D[G1D<1e-7]=0
    G1D=G1D/(max(G1D))*Fr
    Fr.tmp[(1+burns):nt,]=G1D
    
  }
  
  
  return(Fr.tmp)
  
}

#=============================================================================
#Calculate the population spread rate from an IGR. This is based on
#the math from Neubert and Caswell 00, but originating with Weinberger 78, 
#and Kot 92. 
#=============================================================================

get.spread.rate = function(gr1.fast,a_rr,sr) {
  u=seq(0,.1,0.0001)
  #For Lalpacian kernel 
  r1=gr1.fast*1/(1-(1/a_rr[1])^2*u^2) 
  r2=sr[1]*1
  cs=1/u*log(r1+r2)
  cs_all=min(cs[cs>0],na.rm=T)#/10
  
  return(cs_all)
}


#=============================================================================
#Calculate a new fitness distribution for a single time step. 
#=============================================================================

get.fast.peak = function(peak.new,Fr,Ds,cs_all,stc,np){
  peak.new=peak.new+cs_all 
  max.new=Fr[1]
  
  X2.new=stc[(1:(np)),2]
  
  G1D.new=exp(-((X2.new-peak.new)^2)/(2*Ds[1]));
  G1D.new[G1D.new<1e-7]=0
  G1D.new=G1D.new/(max(G1D.new))*max.new
  
  return(G1D.new)
  
}

#=============================================================================
#Population dynamics
#=============================================================================
#One step of the spatial Leslie-Gower model for N species
#1D array of sites. 
#It takes as variables: 
#	Fr.spp 	The intrinsic ranges of all species
#	nr.spp 	Population matrix with current populations or initial conditions
#	sr.spp 	Species survival
#	alpha.spp Competition matrix
#	kd.spp 	Dispersal kernels
#	fkd.yes If this is yes, then it is the Fourier transformed kernels being
#			passed to the function -- Removed for now due to padding!fkd.yes=FALSE,
#	kc.spp 	Competition kernels
#=============================================================================

pop_lg = function (Frs.spp,nr.spp, sr.spp,alpha.spp, kd.spp, kc.spp, pad = FALSE){
  
  
  Frs.spp = as.matrix(Frs.spp)
  nr.spp = as.matrix(nr.spp)
  np_pre = dim(as.matrix(Frs.spp))[1] #Space
  nspp = dim(as.matrix(Frs.spp))[2] #Number of species
  nspp2 = dim(as.matrix(kc.spp))[2]
  
  if (pad == TRUE) {
    padding1 = matrix(0,np_pre,nspp)
    padding2 = matrix(0,np_pre,nspp2)
    Frs.spp = rbind(padding1, Frs.spp,padding1)
    nr.spp = rbind(padding1,nr.spp,padding1)
    kd.spp = rbind(padding2,kd.spp,padding2)
    kc.spp = rbind(padding2,kc.spp,padding2)
    pad1 = np_pre
  } else {pad1 = 0}
  
  np = dim(as.matrix(Frs.spp))[1] #Space
  
  #Check the dimensionality of the dispersal kernel. If it is a spatially 
  #heteroeneous dispersal kernel, then set het_kern=1
  
  het_kern=0
  
  if (!is.na(dim(kd.spp)[3])){ het_kern=1 } 
  
  #When fkd.yes=FALSE: 
  #Transform dispersal and competition kernels
  if( fkd.yes == FALSE){ 
    if(het_kern==1){
      fkd = kd.spp*0
      
      if(dim(kd.spp)[2]<2){ 
        fkc=fft(as.matrix(kc.spp))
        
        for(ss in 1:nspp){
          fkd[,,ss]=fft(as.matrix(kd.spp[,,ss]))
        }
        
      }else{
        fkc=mvfft(as.matrix(kc.spp))
        
        for(ss in 1:nspp){
          fkd[,,ss]=mvfft(as.matrix(kd.spp[,,ss]))
        }
      }
      
    }else{ 
      
      if(dim(as.matrix(kd.spp))[2]<2){ 
        fkd=fft(as.matrix(kd.spp))
        fkc=fft(as.matrix(kc.spp))
        
      }else{
        
        fkd=mvfft(as.matrix(kd.spp))
        fkc=mvfft(as.matrix(kc.spp))
      }
    }
  } else {
    
    fkd = kd.spp
    fkc = kc.spp
  }
  
  
  nr.burns =array(c(matrix(0.0,2,np),matrix(0.0,2,np)),dim=c(2,np,nspp)) 
  nr.burns[1,,] = nr.spp
  #Seedlings germinate and produce seeds at rate Fr, weighted by competition
  if(dim(as.matrix(nr.burns[1,,]))[2]<2){
    fp = as.matrix(fft(nr.burns[1,,]))
  }else{ 
    fp = as.matrix(mvfft(nr.burns[1,,]))
  }
  
  #Total competition for each resident species 
  cr=matrix(0,np,nspp^2)
  cr_tot = matrix(0,np,nspp)
  for (sa in 1:nspp) { 
    for (sb in 1:nspp) {
      
      
      cr[,(sa*(sa-1)+sb)] = Re(fft((alpha.spp[sa,sb]*fp[,sb]*fkc[,sb]),
                                   inverse=T)/(np+1))
      
      #This if/else makes the code more robust to odd/even spatial extent
      if(ceiling(np/2) == floor (np/2) ){
        cr[,(sa*(sa-1)+sb)] = c(cr[((ceiling(np/2)+1):(np)),(sa*(sa-1)+sb)], 
                                cr[(1:floor(np/2)),(sa*(sa-1)+sb)])
      }else {
        cr[,(sa*(sa-1)+sb)] = c(cr[(ceiling(np/2):(np)),(sa*(sa-1)+sb)], 
                                cr[(1:floor(np/2)),(sa*(sa-1)+sb)])
      }
      
      cr_tot[,sa] = cr_tot[,sa] + cr[ ,(sa*(sa-1)+sb)]	
    }
  }
  
  #competition-weighted seed production
  lam = matrix(0,np,nspp)
  
  #Seed dispersal
  fd = matrix(0,np,nspp)
  fd2= array(matrix(0,np,np),dim=c(np,np,nspp) )
  nr_disp = matrix(0,np,nspp)
  
  for(sa in 1:nspp) { 
    #competition-weighted seed production (Leslie-Gower)
    lam[,sa]=(Frs.spp[,sa])/(1+cr_tot[,sa])
    lam[,sa]= lam[,sa]*(nr.burns[1,,sa]>1e-34)
    #If both species have zero germination, then NAs are produced: 
    lam[,sa][is.na(lam[,sa]) ] = 0
    
    #Seeds disperse
    
    if(het_kern==1){
      #The most robust way to do this is to disperse each site individually 
      #according to its numerical dispersal kernel
      lam.tmp = lam[,sa]*nr.burns[1,,sa]
      lam2 = matrix(0,np,np)
      #lam2[ceiling(np/2),] = lam2[ceiling(np/2),]+t(as.matrix(lam.tmp))
      diag(lam2) = as.matrix(lam.tmp)
      fd2[,,sa]=mvfft(lam2)
      nr_disp[,sa] = rowSums( Re(mvfft( (fd2[,,sa]*fkd[,,sa]), inverse=T)/(np+1)),na.rm=T	)
      
    } else { 
      
      fd[,sa]=fft(lam[,sa]*nr.burns[1,,sa])
      nr_disp[,sa] = Re(fft( (fd[,sa]*fkd[,sa]), inverse=T)/(np+1))
    }
    
    #This if/else makes the code more robust to odd/even spatial extent
    if(ceiling(np/2) == floor (np/2) ){
      nr_disp[,sa] = c(nr_disp[((ceiling(np/2)+1):(np)),sa],nr_disp[(1:floor(np/2)),sa])
    }else {
      nr_disp[,sa] = c(nr_disp[(ceiling(np/2):(np)),sa],nr_disp[(1:floor(np/2)),sa])
    }
    
    #Final population step
    nr.burns[2,,sa] = nr_disp[,sa] + nr.burns[1,,sa]*sr[sa]
    
  }
  
  return(nr.burns[2,((pad1+1):(pad1+np_pre)),])
}

#=============================================================================
#Invader and resident equilibrium densities
#=============================================================================

#=======================#
#Multi-species routines
#=======================#

#=============================================================================
#Get the resident equilibrium densities using simulations. This version does
#not require any analytical calculations and works for an arbitrary number 
#of species. 
#1D array of sites. 
#It takes as variables: 
#	Fr.spp 	The intrinsic ranges of all species
#	s.index The identities of the residents
#	sr.spp 	Species survival
#	alpha.spp Competition matrix
#	kd.spp 	Dispersal kernels
#	fkd.yes If this is yes, then it is the Fourier transformed kernels being
#			passed to the function
#	kc.spp 	Competition kernels
#
# 	burns	Number of time steps at which to test equilibrium: default=500
#	burns.max Maximum number of iterations (in multiples of burns)
#			to test for equilibrium. Default = 10
#	tol 	The tolerance for allowing equilibrium, in terms of significance of
#			fit of a line through the last burns/2 points. Default alpha = 0.05
#
#	fast 	If TRUE, the check for equilibrium is switched to a different routine
#			which is much faster, but potentially less accurate. 
#	fast.tol Test for equilibrium: how close is slope to 0? Default = 1e-5
#=============================================================================

#
get.res.eq = function(Frs.spp,s.index,sr.spp,alpha.spp,kd.spp,  kc.spp, pad = FALSE,burns=500, burns.max = 10,tol=1,fast=FALSE,fast.tol=1e-5){
  
  np = dim(as.matrix(Frs.spp))[1] #Space
  nspp = length(s.index) #Number of residents
  
  Frs.spp = Frs.spp[,s.index]
  sr.spp = sr.spp[s.index]
  #alpha.spp = alpha.spp[s.index,s.index]
  
  nr.burns=array(c(matrix(0.0,burns,np),matrix(0.0,burns,np)),dim=c(burns,np,nspp)) 
  #Initialize residents
  for ( sa in 1:nspp) { nr.burns[1,,sa] = matrix(0.01,1,np)}
  b=1
  #Test to see whether the maximum number of iterations has been reached
  while (b <= burns.max){ 
    #Start each burn iteration
    for (n in 1:(burns-1)) {
      #print(paste ("first", n, sep=" "))
      nr.burns[(n+1),,] = pop_lg(Frs.spp ,nr.burns[n,,], sr.spp,alpha.spp,kd.spp,  kc.spp, pad)
      
    }
    
    #Test for equilibrium for each specie: has the equilibrium density stayed unchanged,
    #according to the significance of a linear fit through the last burns/2 points? 
    
    if(fast==FALSE) { 
      is.eq = matrix(0,nspp,1)
      tt=1:(burns/2+1)
      for(sa in 1:nspp) {
        nr.tmp = rowMeans(nr.burns[(burns/2):burns,,sa])
        sp.lm = lm(nr.tmp~tt)
        if( summary(sp.lm)$coefficients[2,4] < tol ) {is.eq[sa] = 1}
      } 
    }else { 
      
      is.eq = matrix(0,nspp,1)
      tt=1:(burns/2+1)
      for(sa in 1:nspp) {
        nr.tmp = rowMeans(nr.burns[(burns/2):burns,,sa])
        sp.lm = lm(nr.tmp~tt)
        if( sp.lm$coefficients[2] > fast.tol ) {is.eq[sa] = 1}
      }
    }
    
    #Test whether all species have reached equilibrium. If yes, quit and return
    #nr.burns. Otherwise, continue. 
    if( sum(is.eq[sa] == 0) ) {b = burns.max+1}else{ 
      nr.burns[1,,] =nr.burns[burns,,]
      b=b+1
    }
    #print(paste ("outer", n, sep=" "))
    
  }
  
  #If this quit due to equilibrium, return equilibriums. Otherwise, return an error
  if (b > burns.max) { return(nr.burns[burns,,])} else {
    return( print("Error: Max iterations exceded without reaching equilibrium"))}
  
}

#============================================================================
#Get the invader's low-density equilibrium density using simulations. This version does
#not require any analytical calculations and works for an arbitrary number 
#of species. 
#1D array of sites. 
#It takes as variables: 
#	Fr.inv	The intrinsic range of the invader
# 	nr.res 	The residents' stationary distributions
#	s.inv 	The identity of the invader
#	sr.inv 	Invader survival
#	alpha.spp Competition matrix
#	kd.spp 	Dispersal kernels
#	fkd.yes If this is yes, then it is the Fourier transformed kernels being
#			passed to the function
#	kc.spp 	Competition kernels
#
# 	burns	Number of time steps at which to test equilibrium: default=500
#	burns.max Maximum number of iterations (in multiples of burns)
#			to test for equilibrium. Default = 10
#	tol 	The tolerance for allowing equilibrium, in terms of significance of
#			fit of a line through the last burns/2 points. Default alpha = 0.05
#
#	fast 	If TRUE, the check for equilibrium is switched to a different routine
#			which is much faster, but potentially less accurate. 
#	fast.tol Test for equilibrium: how close is slope to 0? Default = 1e-5
#=============================================================================

get.inv.ldeq = function(Frs.spp, nr.res,s.inv,sr.spp,alpha.spp,kd.spp,  kc.spp, pad = FALSE,burns=500, burns.max = 10,tol=1,fast=FALSE,fast.tol=1e-5){
  
  np = dim(as.matrix(Frs.spp))[1] #Space
  nspp = dim(as.matrix(Frs.spp))[2] #Number of species
  s.index= 1:nspp
  s.index = s.index[-s.inv]
  
  inv.id =1e-5
  
  nr.burns=array(c(matrix(0.0,burns,np),matrix(0.0,burns,np)),dim=c(burns,np,nspp)) 
  
  #Initialize invader and residents
  nr.burns[1,,s.inv] = matrix(inv.id,1,np)
  nr.burns[1,,s.index] = nr.res[,s.index]
  
  b=1
  #Test to see whether the maximum number of iterations has been reached
  while (b <= burns.max){ 
    #Start each burn iteration
    for (n in 1:(burns-1)) {
      nr.burns[(n+1),,s.inv] = pop_lg(Frs.spp ,nr.burns[n,,], sr.spp,alpha.spp,kd.spp, kc.spp, pad )[,s.inv]
      
      
      #Reset invader to low density
      nr.burns[(n+1),,s.inv] = (inv.id /sum(nr.burns[(n+1),,s.inv])) * nr.burns[(n+1),,s.inv]
      
    }
    
    #Test for equilibrium for invader: has the density stayed unchanged,
    #according to the significance of a linear fit through the last burns/2 points? 
    
    is.eq =0
    tt=1:(burns/2+1)
    if(fast==FALSE) { 
      
      nr.tmp = rowMeans(nr.burns[(burns/2):burns,,s.inv])
      sp.lm = lm(nr.tmp~tt)
      if( summary(sp.lm)$coefficients[2,4] < tol ) {is.eq = 1}
      
    }else { 
      
      nr.tmp = rowMeans(nr.burns[(burns/2):burns,,s.inv])
      sp.lm = lm(nr.tmp~tt)
      if( sp.lm$coefficients[2] > fast.tol ) {is.eq = 1}
      
    }
    
    #Test whether the invader has reached equilibrium. If yes, quit and return
    #nr.burns. Otherwise, continue. 
    if( sum(is.eq == 0) ) {b = burns.max+1}else{ 
      nr.burns[1,,s.inv] =nr.burns[burns,,s.inv]
      b=b+1
    }
    
  }
  
  #If this quit due to equilibrium, return equilibriums. Otherwise, return an error
  if (b > burns.max) { return(nr.burns[burns,,s.inv])} else {
    return( print("Error: Max iterations exceded without reaching equilibrium"))}
  
}
