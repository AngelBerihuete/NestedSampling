NestedSampling <- function(prior, loglikelihood, dim.par = 3, M = 100, N = 100){

  # Functions required
  require(plyr)
    
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
  f <- -2 # As Farhan Feroz in doi:10.1111/j.1365-2966.2007.12353.x (0.2 in log-evidence)
  max.loglik <- max(nested.object$S$loglik)
  logXj <- 0
  cond2 <- (max.loglik + logXj) > f + logZ

  while(cond1 & cond2) {
    evolved.nested.object <- NestedSamplingStep(nested.object, i, prior, loglikelihood)
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

NestedSamplingStep <- function(nested.object, iteration, prior, loglikelihood) {
  
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
  new.active.point <- SamplingNewCandidate(worstloglik, method = 'prior', prior, loglikelihood)
  S[worstpoint,] <- new.active.point
  
  nested.object$S <- S
  nested.object$logZ <- logZnew
  nested.object$H <- H
  
  return(list(no = nested.object, ip = inactive.point))
}


SamplingNewCandidate <- function(worstloglik, method = 'prior', prior, loglikelihood){
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

log.plus <- function(x,y) {
  if(x>y) x+log(1+exp(y-x))
  else    y+log(1+exp(x-y))
}
