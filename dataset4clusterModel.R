# data {
#   // General setup
#   int<lower=1> N; // number of objects observed
#   int<lower=3> D; // dimension of the observations
#   int<lower=0> E;
#   row_vector[D] y_obs[N]; // observations
#   real lambda_obs[N]; // (observarion within the cluster)
#   
#   // Setup the cluster
#   int<lower=0> K; // max order of chebyshev polynomial
#   row_vector[K] coefs_mu_0[E];
#   // Setup the field
# }

require(reshape2)
require(plyr)
require(signal)
require(chebpol)
require(mvtnorm)


# code coefs_mu_0 #needed for stan code
#######################################
load("members-fulldataset.RData")
load('../../pleiades/datasets/iter5.RData')
cm.members <- cm
cov.members <- cov
dens.members <- dens
s.memebers <- s
lambdas.members <- lambdas

load("iter5.RData")

# Normalizo lambda para que esté entre -1 y 1 (aconsejado por el paquete chebpol)
lambdas <- pc$lambda
lambdas <- lambdas - min(lambdas)
lambdas <- lambdas/ max(lambdas)
lambdas <- lambdas/ 0.5
lambdas <- lambdas - 1

# Defino una grid de chebyshev entre -1 y 1 con 8 puntos
K <- 4 # fixed order of chebyshev polynomials
lambdas.out<-chebknots(K)

# Obtengo los valores de la curva principal en esa grid de chebyshev
method="spline"
s<-cbind(
  interp1(lambdas,pc$s[,1],xi=lambdas.out[[1]],method,extrap=TRUE,'natural'),
  interp1(lambdas,pc$s[,2],xi=lambdas.out[[1]],method,extrap=TRUE,'natural'),
  interp1(lambdas,pc$s[,3],xi=lambdas.out[[1]],method,extrap=TRUE,'natural'),
  interp1(lambdas,pc$s[,4],xi=lambdas.out[[1]],method,extrap=TRUE,'natural'),
  interp1(lambdas,pc$s[,5],xi=lambdas.out[[1]],method,extrap=TRUE,'natural')
)

coef <- matrix(NA,nrow=5,4)

# Ahora, para cada variable, descompongo la curva en la base de chebyshev
for (i in 1:5)
{
  k <- chebcoef(s[,i])
  coef[i,] <- k
  tch <- Vectorize(function(x) chebeval(x,k))
}

coefs_mu_0 <- t(coef) # needed for stan code

######################################################

load('../../stage2/final.RData')
load('full-dataset-23.09.2013.RData')




tmp <- sort(lambdas)
lsb.beta <- chebcoef(s[,5])
tch2 <- Vectorize(function(x) chebeval(x,lsb.beta))
estim.y <- tch2(tmp)

load('../../stage2/final.RData')
cnm.no.members <- cnm
cov.no.members <- cov
s.no.members <- s

# rename according with stan code
y_obs <- 
N <- dim(y_obs)[1] # number of observations
D <- dim(y_obs)[2] # dimension of the observations

coefs_mu_0 <- 

dump(c("N", "K","x_obs","y_obs"),"clusterModel.dataset.R")
