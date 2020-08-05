setwd("C:/Users/Rakshanda/Desktop/Bayesian Learning/LABS/lab_03")
data = read.csv("rainfall.dat")
library("LaplacesDemon")
library(latex2exp)
library(coda)
set.seed(12345)
#______(a)_______#
x = data$X136
hist(x, breaks = 50,probability = T)
lines(density(x),col='red',lwd = 1.5, main = 'Distribution of Data')
#___Conditional distributions of mu and sigma2___#
update_mu = function(nDraws,y,sigma2,tau2_0,mu0){
  n = length(y)
  tau2_n = (n/sigma2 + 1/tau2_0)^(-1)
  w = (n/sigma2)/(n/sigma2 + 1/tau2_0)
  mu_n = w*mean(y)+(1-w)*mu0
  mu = rnorm(1,mean = mu_n,sd = sqrt(tau2_n))
  return(mu)
}
###### Defining a function that simulates from the 
rInvchisq <- function(n, df, scale) (df*scale)/rchisq(n,df=df)

update_sigma2 = function(nDraws,y,mu,v0,sig2_0){
  n = length(y)
  vn = v0 + n
  sig2_n = (v0*sig2_0 + sum((y-mu)^2)) / vn
  sigma2  = rInvchisq(n = 1,df = vn,scale = sig2_n)
  return(sigma2)
}
#___Posterior distribution using gibbs___#
gibbs =function(nDraws=1000,y,initial,prior){
  #browser()
  # initial values
  mu = initial[1]
  sigma2 = initial[2]
  
  draws = matrix(nrow = nDraws,ncol = 2)
  draws[1,] = c(mu,sigma2)
  for (i in 2:nDraws) {
    
    # update mu
    mu = update_mu(1,y,sigma2,tau2_0 = prior$tau2_0,mu0 = prior$mu0)
    draws[i,1] = mu
    
    # update sigma2 
    sigma2  = update_sigma2(1,y,mu,v0 = prior$mu0,sig2_0 = prior$sig2_0)
    draws[i,2] = sigma2
  }
  draws
}
#___prior hyperparameters___#
prior = list()
prior$mu0 = 0
prior$tau2_0 = 50
prior$v0 = 10
prior$sig2_0 = 100
n =  length(x)
hist(x,breaks = 50,probability = T, ylim = c(0,0.04))
curve(expr = dnorm(x = x,mean = prior$mu0,sd = sqrt(prior$sig2_0)),
      lty=2,col = 'blue',add = T,lwd=2)
#____________#
init_val = c(mean(x),var(x))
sample = gibbs(nDraws = 2000, y = data$X136,initial = init_val,prior = prior)
#plot(as.mcmc(sample))
summary(as.mcmc(sample))
mu_1 = sample[,1]
sigma2_1 = sample[,2]

#__ Following is the contour plot of the joint posterior of(mu,sigma2)
plot(x = mu_1,y = sigma2_1, xlab = TeX("$\\mu$"), ylab = TeX("$\\sigma^2$"))
abline(h = var(x), v = mean(x), col = 'red')

#____Time series plots____#
ts.plot(mu_1,ylab = TeX(paste("$\\mu$")))

hist(mu_1, xlab = 'Gibbs Sample',probability = T,
     main = TeX(paste("$\\mu |\\sigma^2,X$")),breaks = 30)
lines(density(mu_1), col ='red',lwd = 2)
densplot(as.mcmc(mu_1),  TeX(paste("$\\mu |\\sigma^2,X$")))
#par(mfrow=c(1,2))
ts.plot(sigma2_1,ylab = TeX(paste("$\\sigma^2$")))
hist(sigma2_1, xlab = 'Gibbs Sample',probability = T,
     main = TeX(paste("$\\sigma^2 |\\mu,X$")), breaks = 30)
lines(density(sigma2_1), col ='red',lwd = 2)

#____Convergence using gelman rubin diagnostic____#
set.seed(12345)
init= c(0,1)
post1 = gibbs(nDraws = 2000,y = x,initial = init,prior)
init= c(5,0.01)
post2 = gibbs(nDraws = 2000,y = x,initial = init,prior)
init= c(10,0.5)
post3 = gibbs(nDraws = 2000,y = x,initial = init,prior)
init= c(20,0.4)
post4 = gibbs(nDraws = 2000,y = x,initial = init,prior)
init= c(15,10)
post5 = gibbs(nDraws = 2000,y = x,initial = init,prior)

pmc_mu = mcmc.list(as.mcmc(post1[,1]),as.mcmc(post2[,1]),as.mcmc(post3[,1]),
                   as.mcmc(post4[,1]),as.mcmc(post5[,1]))

pmc_sig = mcmc.list(as.mcmc(post1[,2]),as.mcmc(post2[,2]),as.mcmc(post3[,2]),
                   as.mcmc(post4[,2]),as.mcmc(post5[,2]))

a = gelman.diag(x = pmc_mu)
b = gelman.diag(x = pmc_sig)
#par(mfrow=c(1,2))
gelman.plot(pmc_mu)
gelman.plot(pmc_sig)
 
#______(b)______#
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

xGrid = seq(min(x)-1*apply(data,2,sd),max(x)+1*apply(data,2,sd),length = 100)

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
  
  #Computing mixture density means or thetas
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

t = gibbs_mixture(x,nDraws = 1000,K = 2,prior = prior)
mu = t$mu
sigma2 = t$sigma2
theta = t$theta
# conditional distributions
par(mfrow=c(1,2))
hist(mu[,1], xlab = 'Gibbs Sample',probability = T,ylim = c(0,1.2),xlim = c(7,15),
     main = TeX(paste("$\\mu_1|I,\\sigma^2_1,X$")),breaks = 100)
lines(density(mu[,1]), col ='red',lwd = 2)

hist(mu[,2], xlab = 'Gibbs Sample',probability = T,xlim = c(50,65),
     main = TeX(paste("$\\mu_2|I,\\sigma^2_2,X$")),breaks = 100)
lines(density(mu[,2]), col ='red',lwd = 2)
#______________________#
hist(sigma2[,1], xlab = 'Gibbs Sample',probability = T,
     xlim = c(30,150),ylim = c(0,0.07),
     main = TeX(paste("$\\sigma^2_1|I,\\mu_1,X$")),breaks = 150)
lines(density(sigma2[,1]), col ='red',lwd = 2)

hist(sigma2[,2], xlab = 'Gibbs Sample',probability = T, xlim = c(1600,2500),
     main = TeX(paste("$\\sigma^2_2|I,\\mu_2,X$")),breaks = 50)
lines(density(sigma2[,2]), col ='red',lwd = 2)
#___convergence plots___#
library(latex2exp)
plot(mu[,1], type = 'l', ylim = c(0,70),ylab = TeX("$\\mu$"))
lines(mu[,2], col = 'red')
legend("topright", box.lty = 1, bty = 'n',
       legend = TeX(c("$\\mu_1$","$\\mu_2$")), col=c("black","red"), lwd = 2)

plot(sigma2[,1], type = 'l',ylim = c(0,3000), ylab = TeX("$\\sigma^2$"))
lines(sigma2[,2], col = 'red')
legend("topright", box.lty = 1, bty = 'n', col=c("black","red"), lwd = 2,
       legend = TeX(c("$\\sigma_1^2$","$\\sigma_2^2$")))
#_____________#
hist(x, breaks = 50, probability = T,
     xlim = c(min(xGrid),max(xGrid)),
     main = "Final Fitted Density")
lines(x = xGrid, y = theta, col = 'red', lwd = 2)
lines(xGrid,dnorm(xGrid, mean = mean(mu_1[c(1000:2000)]),sd = sqrt(mean(sigma2_1[c(1000:2000)]))), 
      col = 'blue', lwd = 2, lty = 2)
legend('topright', legend = c("Guassian Mixture Model","Normal Model"),
       col=c("red","blue"), lwd = 1, bty = 'n', lty = c(1,2))


