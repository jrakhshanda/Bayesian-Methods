setwd("C:/Users/Rakshanda/Desktop/Bayesian Learning/LABS/lab_03")
data = read.csv('rainfall.dat')
x = as.matrix(data$X136)

set.seed(12345)

rInvchisq = function(n, df, scale) (df*scale)/rchisq(n,df=df)

####### Defining a function that simulates from a Dirichlet distribution for pi
rDirichlet = function(alpha){
  n = length(alpha)
  x = rep(0,n)
  for (i in 1:n){
    x[i] = rgamma(1,alpha[i],1)
  }
  x = x/sum(x)
  return(x)
}

# Simple function that converts between two different representations of the mixture 
# allocation. 
Alloc = function(I){
  n = nrow(I)
  vect = rep(0,n)
  for (i in 1:n){
    vect[i] = which(I[i,] == 1)
  }
  return(vect)
}

#' @param nDraws number of samples
#' @param x given data
#' @param K number of components for mixture model
#' @param prior a list of hyperparameters of priors
#' @param prior$alpha Dirichlet(alpha)
#' @param prior$mu0 prior mean of mu
#' @param prior$tau2_0 prior SD of mu
#' @param prior$sigma2_0 prior or best guess of sigma 2
#' @param prior$v0 degrees of freedom for prior on sigma2


gibbs_mixture = function(nDraws=100,x,K=2,prior)
{
  n = length(x)
  
  # Initial mcmc values
  I = t(rmultinom(n, size = 1 , prob = rep(1/K,K))) 
  mu = rep(0,K)
  sigma2 = rep(var(x),K)
  probObsInComp = rep(NA,K)
  
  #variables to store final values
  MU = matrix(nrow = nDraws,ncol = K)
  Sigma2 = matrix(nrow = nDraws,ncol = K)
  
  #Computing mixture density means or thetas and plotting
  xGrid = seq(min(x)-1*apply(data,2,sd),max(x)+1*apply(data,2,sd),length = 100)
   mixDensMean = rep(0,length(xGrid))
   count = 0
  
  for (it in 1:nDraws)
  {
    alloc = Alloc(I)
    nAlloc = colSums(I)
    
    # Update components probabilities
    pi = rDirichlet(prior$alpha + nAlloc)
   
     # Update mu's
    for (j in 1:K){
      precPrior = 1/prior$tau2_0[j]
      precData = nAlloc[j]/sigma2[j]
      precPost = precPrior + precData
      w = precPrior/precPost
      mu_n = w*prior$mu0 + (1-w)*mean(x[alloc == j])
      tau2_n = 1/precPost
      mu[j] = rnorm(1, mean = mu_n, sd = sqrt(tau2_n))
    }
    MU[it,] = mu 
   
    # Update sigma2's
    for (j in 1:K){
      sigma2[j] = rInvchisq(1, df = prior$v0[j] + nAlloc[j], 
                scale = (prior$v0[j]*prior$sigma2_0[j] + 
                        sum((x[alloc == j] - mu[j])^2)) / (prior$v0[j] + nAlloc[j]))
    }
    Sigma2[it,] = sigma2
    
    # Update probability Observations in Computation
    for (i in 1:n){
      for (j in 1:K){
        probObsInComp[j] = pi[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
      }
      I[i,] = t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
    }
    
    # mixture of densities 
    if (TRUE && (it%%1 ==0)){
      mixDens = rep(0,length(xGrid))
      count = count + 1
      for (j in 1:K){
        temp = pi[j] * dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
        mixDens = mixDens + temp
      }
      mixDensMean = ((count-1)*mixDensMean + mixDens)/count
    }
  }
   
  lst = list(mu = MU, sigma2 = Sigma2, theta = mixDensMean)
  return(lst)
}
# Prior options
prior = list()
prior$alpha = 10*rep(1,2) 
prior$mu0 = rep(0,2) 
prior$tau2_0 = rep(10,2)
prior$sigma2_0 = rep(var(x),2) 
prior$v0 = rep(4,2) 

t = gibbs_mixture(x = data$X136,nDraws = 100,K = 2,prior = prior)

mu = t$mu
sigma2 = t$sigma2
theta = t$theta 

#___convergence plots___#
library(latex2exp)
plot(mu[,1], type = 'l', ylim = c(0,70))
lines(mu[,2], col = 'red')
legend("topright", box.lty = 1, bty = 'n',
       legend = TeX(c("$\\mu_1$","$\\mu_2$")), col=c("black","red"), lwd = 2)

plot(sigma2[,1], type = 'l',ylim = c(0,2500))
lines(sigma2[,2], col = 'red')
legend("topright", box.lty = 1, bty = 'n', col=c("black","red"), lwd = 2,
       legend = TeX(c("$\\sigma_1^2$","$\\sigma_2^2$")))
#_____________#
hist(x, breaks = 50, freq = F,probability = T,
     xlim = c(min(xGrid),max(xGrid)), main ='Histogram of Data')
lines(x = xGrid, y = theta, col = 'red', lwd = 2, lty = 2)

