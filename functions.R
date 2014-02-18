# Dataset: uniform in angle and radio uniform entre 3 y 4.

NestedSampling <- function(prior, loglikelihood, dim.par = 3, M = 100, N = 100){
  # Functions required
  require(plyr)
  log.plus <- function(x,y) {
    if(x>y) x+log(1+exp(y-x))
    else    y+log(1+exp(x-y))
  }
  
  # Creating the first bunch of alive points
  S <- matrix(NA, N, dim.par)
  for(j in 1:N){
    S[j,]  <- prior()
  }
  loglik <- adply(S, .margins = 1, loglikelihood)
  S <- as.data.frame(S)
  S$loglik <- loglik$V1
  
  # Setting up the logZ and H, then nested.object
  H    <- 0 # Information, initially 0
  logZ <- -.Machine$double.xmax # ln(Evidence Z, initially 0)
  nested.object <- list(S = S, logZ = logZ, H = H)
  inactive.points <- data.frame()
  # setting up loop termination
  i <- 1
  cond1 <- i < M
  f <- -2 # As Farhan Feroz in doi:10.1111/j.1365-2966.2007.12353.x
  max.loglik <- max(nested.object$S$loglik)
  logXj <- 0
  cond2 <- (max.loglik + logXj) > f + logZ
  
  while(cond1 & cond2) {
    evolved.nested.object <- NestedSamplingStep(nested.object, i, prior, likelihood)
    nested.object <- evolved.nested.object$no
    inactive.points <- rbind(inactive.points, evolved.nested.object$ip)
    # setting up termination conditions
    
    max.loglik <- max(nested.object$S$loglik)
    logXj <- -i/N
    logZ <- nested.object$logZ
    cond2 <- (max.loglik + logXj) > f + logZ
    
    i <- i + 1
    cond1 <- i < M
  }
  return(list(no = nested.object, ip=inactive.points))
  
}

NestedSamplingStep <- function(nested.object, iteration, prior, likelihood) {
  
  S <- nested.object$S
  H <- nested.object$H
  logZ <- nested.object$logZ
  
  worstloglik <- min(S$loglik)
  worstpoint <- which.min(S$loglik)  
  
  inactive.point <- S[worstpoint,]
  
  # Updating nested object
  N <- dim(S)[1]
  logwidth <- log(0.5*(exp(-(iteration-1)/N)-exp(-(iteration+1)/N)))
  logWt <- logwidth + worstloglik
  logZnew <- log.plus(logZ, logWt)
  H <- exp(logWt - logZnew) * worstloglik + exp(logZ - logZnew) * (H + logZ) - logZnew
  
  new.active.point <- SamplingNewCandidate(worstloglik, method = 'prior', prior, likelihood)
  S[worstpoint,] <- new.active.point
  
  nested.object$S <- S
  nested.object$logZ <- logZnew
  nested.object$H <- H
  
  return(list(no = nested.object, ip = inactive.point))
}


SamplingNewCandidate <- function(worstloglik, method = 'prior', prior, likelihood){
  if(method == 'prior'){
    theta <- prior()
    loglik <- loglikelihood(theta)
    
    if(loglik > worstloglik){
      return(c(theta, loglik))
    }else{
      cond <- TRUE
      while(cond){
        theta <- prior()
        loglik <- loglikelihood(theta)
        cond <- loglik < worstloglik
      }
      return(c(theta, loglik))
    }
  }
}

prior <- function(){
  u <- runif(3)
  theta <- 20*u-10
  return(theta)
}

loglikelihood <- function(theta) -0.5*sum(theta^2)

prior2 <- function() qunif(runif(2),min=0,max=10*pi)
loglikelihood2 <- function(theta)  (2+cos(theta[1]/2)*cos(theta[2]/2))^5



DoMeasurements <- function(results){
  N <- dim(results$no$S)
  logZ=results$no$logZ
  logZ.sd=sqrt(results$no$H/N)
  Hnats=results$no$H
  Hbits=results$no$H/log(2)
  return(c(logZ, logZ.sd, Hnats, Hbits))
}

DoMeasurements(results)
ggplot(data = results$ip, aes(V3)) + geom_density()
ggplot(data = results$ip, aes(V2,V3)) + geom_density2d() + geom_point(aes(colour = loglik))

# -----------------------------------
# Optimal decomposition for Multinest
# -----------------------------------
