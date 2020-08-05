setwd("C:/Users/Rakshanda/Desktop/Bayesian Learning/LABS/lab_03")
data = read.table("eBay.dat",header = T)
library(mvtnorm)
library(coda)
options("scipen"=100)
set.seed(123456)
y = data$nBids
X = as.matrix(data[,-1])
#_______(a)______#
model = glm(formula = nBids ~0+.,data = data,family = "poisson")
summary(model)
MLE_coeff = matrix(round(model$coefficients,),nrow = 1)
colnames(MLE_coeff) = colnames(X)

#_____(b)______#
#' @param beta regression coefficients
#' @param y response variable
#' @param X features/independen variables
#' @param sigma prior or best guess of sigma 2

# prior setting
Sigma = 100 * solve(t(X)%*%X)
mu0 = rep(0,ncol(X))

log_posterior = function(beta,y,X,Sigma,mu0){
  beta = as.vector(beta)
  linPred = X%*%beta
  log_lik = sum(y * linPred) - sum(exp(linPred))
  if (abs(log_lik) == Inf) log_lik = -20000
  log_prior = dmvnorm(beta, matrix(c(mu0),length(beta),1), Sigma, log=TRUE)
  return(log_lik + log_prior)
}
optimal = optim(rep(0,9),log_posterior,gr=NULL,y,X,Sigma,mu0,method=c("BFGS"),
                control=list(fnscale=-1),hessian=TRUE)
betas = optimal$par
jInv = solve(-1*optimal$hessian)

beta = matrix(round(betas,4),nrow = 1)
colnames(beta) = colnames(X)
#______(c)_______#
#' @param postFun Log posterior of the objective distribution
#' @param nDraws random Sample to draw from metropolis hasting sample 
#' @param y response variable
#' @param X features/independen variables
#' @param theta_0 prior or best guess of the parameter
#' @param c tuning parameter
MH = function(postFun, theta_0, c,iter,warmup,...)
{
  nDraws = iter + warmup
  optimal = optim(theta_0, postFun, gr=NULL,...,
                  method=c("BFGS"), control=list(fnscale=-1),
                  hessian=TRUE)
  jInv = solve(-optimal$hessian)*c
  X = matrix(nrow = nDraws, ncol = length(theta_0)) 
  X[1,] = theta_0
  accepted = 0
  count = 0
  for (i in 2:nDraws)
  {
    while(accepted < nDraws) {
     current = X[i-1,]
     theta_p = rmvnorm(1, mean = current, sigma = jInv)
    
    # acceptance probability ratio
    ratio = exp(postFun(c(theta_p),...) - postFun(c(current),...)) *
      dmvnorm(current,theta_p,jInv) / dmvnorm(theta_p,current,jInv)
    
    alpha =  min(1, ratio)
    u = runif(1)
    count = count + 1
    if(u<=alpha){ 
      X[i,] = theta_p
      current = theta_p
      accepted = accepted+1
      break
      }
    }
  }
  return(list(draws = X[-c(1:warmup),], acceptance_rate = accepted/count))
}
smp = MH(log_posterior,rep(0,9), c=0.6,iter = 1000, warmup = 500, y,X, Sigma,mu0)
betas = smp[[1]]

traceplot(as.mcmc(betas))
# windows()
# plot(as.mcmc(betas))
#______(d)______#
#set.seed(123456)
auction = matrix(c(1,1,1,1,0,0,0,1,0.5), nrow = 1)
colnames(auction) = colnames(X)

predictions = function(nDraws=1000, X = auction){
    draws = c()
    for (i in 1:nDraws) {
      probs = exp(X %*% betas[i,])
      draws[i] = rpois(n = 1,lambda = probs)
    }
    draws
}
pred = predictions()
table(pred)
hist(pred, breaks = 40, probability = T, col = "dodgerblue4",
     xlab = 'Predictions', main = 'Histogram of Predictions')

