source('functions.R')

prior <- function(){
  u <- runif(1)
  theta <- 20*u-10
  return(theta)
}

loglikelihood <- function(theta) -0.5*sum(theta^2)

# Generate data

results <- NestedSampling(prior = prior, loglikelihood = loglikelihood, 
                          dim.par = 1, M = 100, N = 100)


prior2 <- function() qunif(runif(2),min=0,max=10*pi)
loglikelihood2 <- function(theta)  (2+cos(theta[1]/2)*cos(theta[2]/2))^5

library(reshape2)
library(ggplot2)
pp <- function (n,r=10) {
  x <- seq(0, r*pi, len=n)
  df <- expand.grid(x=x, y=x)
  df$z <- (2+cos(df$x/2)*cos(df$y/2))^5
  df
}

dataset <- pp(100)
ggplot(data = dataset, aes(x,y)) + geom_tile(aes(fill = z)) 

results <- NestedSampling(prior = prior2, loglikelihood = loglikelihood2,
                          dim.par = 2, M = 1000, N = 100)

ggplot(data = results$no$S, aes(V1,V2)) + geom_point(aes(color = loglik))
ggplot(data = results$ip, aes(V1,V2)) + geom_point(aes(color = loglik))

require(signal)
require(chebpol)
require(mvtnorm)


# Voy a utilizar los datos de siempre para obtener la descomposición de la curva
# principal en la base de polinomios de chebyshev
load("iter5.RData")

# Normalizo lambda para que esté entre -1 y 1 (aconsejado por el paquete chebpol)
lambdas <- pc$lambda
lambdas <- lambdas - min(lambdas)
lambdas <- lambdas/ max(lambdas)
lambdas <- lambdas/ 0.5
lambdas <- lambdas - 1

# Defino una grid de chebyshev entre -1 y 1 con 8 puntos
lambdas.out<-chebknots(8)

# Obtengo los valores de la curva principal en esa grid de chebyshev
method="spline"
s<-cbind(
  interp1(lambdas,pc$s[,1],xi=lambdas.out[[1]],method,extrap=TRUE,'natural'),
  interp1(lambdas,pc$s[,2],xi=lambdas.out[[1]],method,extrap=TRUE,'natural'),
  interp1(lambdas,pc$s[,3],xi=lambdas.out[[1]],method,extrap=TRUE,'natural'),
  interp1(lambdas,pc$s[,4],xi=lambdas.out[[1]],method,extrap=TRUE,'natural'),
  interp1(lambdas,pc$s[,5],xi=lambdas.out[[1]],method,extrap=TRUE,'natural')
)

coef <- matrix(NA,nrow=5,8)

# Ahora, para cada variable, descompongo la curva en la base de chebyshev
for (i in 1:5)
{
  k <- chebcoef(s[,i])
  coef[i,] <- k
  tch <- Vectorize(function(x) chebeval(x,k))
  
  # Por si quieres representar el ajuste 
  # plot(lambdas,pc$s[,i],ty="p",pch=16,cex=.3)
  # tmp <- sort(lambdas)
  # lines(tmp,col="orange",lw=2)
}

# observed lambdas and estimated y with chebyshev polynomials

tmp <- sort(lambdas)
lsb.beta <- chebcoef(s[,1])
tch2 <- Vectorize(function(x) chebeval(x,lsb.beta))
estim.y <- tch2(tmp)


prior3 <- function() {
  theta <- array(rmvnorm(1, mean = coef[1,]))
  return(theta)
}

loglikelihood3 <- function(theta){
  theta <- array(theta)
  tch1 <- Vectorize(function(x) chebeval(x,theta))
  true.y <- tch1(tmp)
  loglik <- mahalanobis(true.y, center = estim.y, cov = diag(length(estim.y)))
  return(loglik)
}

results <- NestedSampling(prior = prior3, loglikelihood = loglikelihood3,
                          dim.par = 8, M = 100, N = 100)


ggplot(data = results$no$S, aes(V1)) + geom_histogram()
ggplot(data = results$no$S, aes(V2)) + geom_histogram()
ggplot(data = results$no$S, aes(V3)) + geom_histogram()




