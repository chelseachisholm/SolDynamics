post.plot <- function(posterior, param=NULL, 
                      param.n=NULL, order.t = F,
                      xmax=NULL,
                      xmin=NULL, 
                      ymin=1, ymax=NULL, ab=0, 
                      y.lab=T, p.dens=T, id.lines=F,
                      dev = F,
                      add.post = NULL, ...){
  #   Posterior is the matrix of the posterior (use as.matrix() if it is an mcmc object)
  #   Param.to.plot is used to select the parameters to plot, defaults to all. 
  #   param.n is the names used 
  #   xmax, xmin, ymin and max are used for plotting, min and max of axes
  #   ab = number or F, if number it is the position of the abline. 
  #   y.lab - do you want the y axis to have labels
  #   p.dens = T/F - Plot densities of posteriors. They are all scaled to prevent overlap.
  #   You can run it by using source_url from the devtools library.
  #	You can post additional posteriors by including them as a list in add.post command
  #   source_url("https://github.com/reuning/EVDebs/blob/master/Bayesian/Posterior-Plots.R")
  
  ##### get defaults and convert to useable things ####
  
  convert.post <- function(tmp.post){
    post.cl <- attributes(tmp.post)$class
    
    if(typeof(tmp.post)=="list"){
      if(post.cl=="rjags") {
        tmp.post <- as.matrix(as.mcmc(tmp.post))
      } else if(post.cl=="mcmc.list") {
        tmp.post <- as.matrix(tmp.post)
      }  else stop("Posterior type not recognized, use either a matrix, jags object, mcmc.list or stanfit")
    } else if(post.cl[[1]]=="stanfit"){
      tmp.post <- as.matrix(tmp.post)
    }
    return(tmp.post)
  }
  
  plot.ci <- function(tmp, p.dens.tmp, y.adj=0, id.lines.tmp=F, type=1){
    ### Function used to plot CI 
    points(y=ii+y.adj, x=mean(tmp))
    ci <- quantile(tmp, c(0.025, .975))
    lines(y=c(ii+y.adj,ii+y.adj), x=ci, lwd=2, lty=type)
    if(id.lines.tmp){
      lines(y=c(ii+y.adj,ii+y.adj), x=c(xmin,xmax), lty=3, lwd=.5)
    }
    if(p.dens.tmp){
      den <- density(tmp)
      lines(y=den$y*conv+ii+y.adj, x=den$x, lty=2) 
    }
  }
  posterior <- convert.post(posterior)
  
  if(is.null(param)){
    param <- colnames(posterior)
    if(!dev & ("deviance" %in% param)){
      param <- param[-which("deviance"==param)]
    }
  } else {
    tmp.param <- character()
    for(ii in 1:length(param)){
      tmp <- grep(param[ii], colnames(posterior), value=T, fixed=T)
      tmp.param <- append(tmp, tmp.param)
      rm(tmp)
    }
    param <- rev(unique(tmp.param))
  }	
  if(is.null(param.n)){
    param.n <- param
  }
  if(ncol(posterior)<length(param)) stop ("More parameters to plot than in posterior")
  if(length(param)!=length(param.n)) stop ("Length of param names not equal to length of params")  
  
  
  add.n <- 0
  if(!is.null(add.post)){
    if(typeof(add.post)=="list"|typeof(add.post)=="S4"){
      if(is.null(attr(add.post, "class"))){
        #### this is a clean list ###
        add.n <- length(add.post)
        
        if(add.n == 1){
          pr.post1 <- convert.post(add.post)
        } else if(add.n == 2){
          pr.post1 <- convert.post(add.post[[1]])
          pr.post2 <- convert.post(add.post[[2]])
        } 
        if(add.n > 2) stop ("Only two additional posteriors can be added") 
      } else {
        pr.post1 <- convert.post(add.post)
        add.n <- 1
      }
    }
  }
  
  
  if(is.null(xmax)) xmax <- max(posterior[,param])
  if(is.null(xmin)) xmin <- min(posterior[,param])
  y.adj <- .1 + .1*add.n
  if(is.null(ymax)) ymax <- length(param)+.1
  
  par(...)
  plot.new()
  plot.window(xlim=c(xmin, xmax),ylim=c(ymin, ymax))
  if(add.n == 0){
    conv <- 1/(1.05*max(apply(posterior[,param], 2, function(ii) max(density(ii)$y))))
  } else if(add.n==1){
    conv <- .9/(1.05*max(apply(posterior[,param], 2, function(ii) max(density(ii)$y))))
  } else {
    conv <- .8/(1.05*max(apply(posterior[,param], 2, function(ii) max(density(ii)$y))))
  }
  
  if(order.t){
    od <- order(apply(posterior[,param], 2, mean))
  } else { od <- length(param):1 }
  for(ii in length(param):1){
    tmp <- posterior[,param[od[ii]]]
    plot.ci(tmp, p.dens.tmp=p.dens, id.lines.tmp=id.lines)
    if(add.n>0){
      tmp.add <- pr.post1[,param[od[ii]]]
      plot.ci(tmp.add, p.dens.tmp=F, y.adj=-.15, id.lines.tmp=F, type=2)
    }
    if(add.n>1) {
      tmp.add <- pr.post2[,param[od[ii]]]
      plot.ci(tmp.add, p.dens.tmp=F, y.adj=-.3, id.lines.tmp=F, type=4)
    }
  }
  xlimlo <- floor(xmin)
  xlimup <- ceiling(xmax)
  n.ax <- length(xlimlo:xlimup)
  if(n.ax < 5){
    x.ax <- xlimlo:xlimup
  } else {
    n.by <- n.ax %/% 5 
    x.ax <- seq(xlimlo, xlimup, by = n.by)
    
  }
  axis(1, at=x.ax)
  if(y.lab){
    axis(2, at=1:length(param), param.n[od], las=1)
  } else {
    axis(2, at=1:length(param), label=NA, tck=0)
  }
  if(ab | ab==0){
    abline(v=ab)
  }
}