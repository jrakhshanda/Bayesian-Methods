options("scipen"=100, "digits"=8)
setwd("C:/Users/Rakshanda/Desktop/Bayesian Learning/LABS/lab_02")
library(mvtnorm)
library(readr)
library(latex2exp)
library(ggplot2)
library(Matrix)
library(matrixcalc)
library(shape)
data = read.table(file = "WomenWork.dat", header = T)

#__________(a)________#
glmModel = glm(Work ~ 0 + ., data = data, family = binomial)
sm = summary(glmModel)
fit = predict(glmModel,newdata = data, type = 'response')
y_hat = ifelse(fit > 0.5, 1L, 0L)
conf_mat = table(data$Work, y_hat)
n = nrow(data)
error_rate = (1-sum(diag(conf_mat))/n)*100
#________(b)__________#
y = as.vector(data$Work)
X = as.matrix(data[,2:9])
nPara = ncol(X)
#___Setting up Log-prior___#
tau = 10
Sigma = tau^2*diag(nPara)

#___Log Posterior Logistic____#
logistic_postrior = function(beta,y,X,sigma){
  nPara = length(beta)
  linPred = X%*%beta
  log_lik = sum( linPred*y -log(1 + exp(linPred)))
  if (abs(log_lik) == Inf) log_lik = -20000
  log_prior = dmvnorm(beta, matrix(0,nPara,1), Sigma, log=TRUE)
  return(log_lik + log_prior)
}
Optimal = optim(rep(0,8),logistic_postrior,gr=NULL,y,X,sigma,method=c("BFGS"),
                   control=list(fnscale=-1),hessian=TRUE)
options("scipen"=100, "digits"=8)
mu = Optimal$par
Sigma = solve(-1*Optimal$hessian)
#___marginal posterior distribution of each beta____#

par(mfrow = c(2,2))
for (i in 2:ncol(X)) {
  grid = seq(mu[i]-3*Sigma[i,i],mu[i]+3*Sigma[i,i],length.out = 1000)
  plot(grid,dnorm(index, mean = mu[i], sd = sqrt(Sigma[i,i])),
       lwd = 2, type = 'l', ylab = 'Density',main = paste("beta",i-1))
}
#_____(a)______#
#____Aprroximation of 95% CI for NSmallChild____#
set.seed(12345678)
beta_7 = rnorm(n = 10000,mean = mu[7],sd = sqrt(Sigma[7,7]))
ci = quantile(beta_7,probs = c(0.025,0.975))
#________(C)________#
options("scipen"=100, "digits"=5)
work_data = matrix(c(1,10,8,10,(10/10)^2, 40, 1,1), nrow = 1)
colnames(work_data) = colnames(X)
predictions = function(size=1000, data = work_data){
  draws = c()
  for (i in 1:size) {
    beta = rmvnorm(n = 1, mean = mu,sigma = Sigma)
    linePred = data %*% t(beta)
    probs = exp(linePred)/(1+exp(linePred))
    draws[i] = rbinom(n = 1,size = 1,probs)
  }
    draws
}
r = predictions(data = work_data)
table(r)
hist(r, main='', xlab = 'Predictions', probability = T)


