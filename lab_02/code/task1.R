setwd("C:/Users/Rakshanda/Desktop/Bayesian Learning/LABS/lab_02")
library(mvtnorm)
library(readr)
library(latex2exp)
library(ggplot2)
library(shape)
data = read.table(file = "TempLinkoping.txt", header = T)
#_______(a)_____#
y = as.matrix(data$temp)
X = matrix(c(data$time^0, data$time, data$time^2), byrow = FALSE, ncol = 3)

#fitting a model with ordinary least squares
model = lm(temp~poly(time,7),data)
sm = summary(model)
ggplot(data, aes(time,temp))+geom_point()+theme_classic()+
  geom_smooth(method = 'lm', formula = y~poly(x,7))
# histogram of residuals
hist(sm$residuals, xlab = 'Residuals', breaks = 10, probability = T)
lines(density(sm$residuals), col = 'red', lwd = 2)
#____conjugate prior____#
rInvchisq <- function(n, df, scale) (df*scale)/rchisq(n,df=df)
n = nrow(data)
mu = c(-10,100,-100)
prior = function(v0 = 4, mu0, sigma20 = 1, omega0= 0.01, size=100)
{
  n = nrow(data)
  reg_curve = matrix(nrow = n, ncol = size)
  for (i in 1:size) {
    omega0 = diag(3)*omega0
    sigma2 = rInvchisq(n = 1,df = v0,scale = sigma20)
    beta = rmvnorm(1, mean = mu0, sigma = sigma2*solve(omega0))
    reg_curve[,i] = X %*% t(beta)
  }
  plot(data$time, y, xlab="Time", ylab="Temperature in Celcius",
       main = TeX(paste0("$\\Omega_0 \ = \ $",omega0, "$\\I_3, \\v_0 \ = \ $",v0,
                         "$, \\sigma_0^2 \ = \ $",sigma20)))
  for (i in 1:size) {
    lines(data$time,reg_curve[,i],type = "l", col=rgb(0.9,0,0,0.2))
  }
}
prior(mu0 = mu)
prior(mu0 = mu, omega0 = 1, v0 = 5,sigma20 = 5)
prior(mu0 = mu, omega0 = 5, v0 = 10,sigma20 = 10)
prior(mu0 = mu, omega0 = 1, v0 = 15,sigma20 = 0.5)
prior(mu0 = mu, omega0 = 5, v0 = 20,sigma20 = 30)

#________(b)______#
mu0 = matrix(c(-10,100,-100), nrow = 3, ncol = 1)
omega0 = 1
v0 = 10
sigma20 = 5

#_____Joint posterior distribution______#
BayesLinReg = function(Y, X, mu0, omega0, v0, sigma20, nIter){
  n = nrow(X)
  k = ncol(X)
  omega0 = diag(k)*omega0
  beta_hat = solve(t(X)%*%X) %*% t(X)%*%Y
  mu_n = solve(t(X)%*%X + omega0) %*% (t(X)%*%X %*% beta_hat+omega0%*%mu0)
  vn = v0+n
  omega_n = t(X)%*%X + omega0
  sigma2_n = (v0*sigma20)+ t(Y)%*%Y + t(mu0)%*%omega0%*%mu0 - t(mu_n)%*%omega_n%*%mu_n
  
  variance = c()
  betas = matrix(nrow = nIter, ncol = ncol(X))
  for (i in 1:nIter) {
    sigma2 = rInvchisq(n = 1,df = vn,scale = sigma2_n/vn)
    variance[i] = sigma2
    beta = rmvnorm(1, mean = mu_n, sigma = c(sigma2)*solve(omega_n))
    betas[i,] = beta
  }
  colnames(betas) = colnames(X)
  return(list(betas=betas, variance=variance))
}
posterior = BayesLinReg(Y = y, X, mu0,omega0, v0, sigma20, nIter = 1000)
betas = posterior[[1]]
sigma2 = posterior[[2]]
#____Marginal Distributions______#
hist(sigma2, breaks = 30,probability = T,xlab = 'Sample',
     main = TeX(paste0("Histogram of $\\sigma^2|y $")))
lines(density(sigma2), col = 'red', lwd = 2)
#____________#
par(mfrow=c(1,3))
for (i in 1:ncol(betas)) {
  hist(betas[,i], breaks = 30,probability = T,xlab = 'Sample',
       main = TeX(paste("Histogram of $\\beta_$",i,"$|y $")))
  lines(density(betas[,i]), col = 'red', lwd = 2)
}

regCurves = matrix(0, nrow = nrow(X), ncol = nrow(betas)) 
for (i in 1:nrow(regCurves)) {
  
  regCurves[i,] = betas%*%c(X[i,]) + sqrt(sigma2)*rnorm(1000,0,1)
}
Median = apply(regCurves,1,median)
lower = apply(regCurves,1,quantile,0.025)
upper= apply(regCurves,1,quantile,0.975)

df = cbind(data, Median,lower,upper)

ggplot(df)+ geom_point(mapping = aes(x = time, y = temp))+theme_classic()+
  geom_line(aes(x = time, y = Median, col = 'Median'), size = 1.5)+
  geom_ribbon(mapping = aes(x= time, ymin= lower, ymax =upper),alpha = 0.4)

#____(c. Highest Temperature)_____#
maxima = mean(-posterior[,2]/(2*posterior[,3]))

plot(data$time, Y, xlab = 'time', ylab = 'temp')
abline(v = maxima, col = 'red')

#_________(d)_________#
x = matrix(c(data$time^0, data$time, data$time^2,data$time^3,
             data$time^4,data$time^5, data$time^6,data$time^7),
             byrow = FALSE, ncol = 8)
mu_7 = c(-12,108, -85,-36, 26,-43,-25,74)
prior_7 = function(mu = mu_7,v0 = 358, sigma20 = 1, lambda= 7, size=100)
{
  n = nrow(data)
  reg_curve = matrix(nrow = n, ncol = size)
  for (i in 1:size) {
    chi = rchisq(n = 1, df = v0)
    tau2 = sum((log(data$time)-mean(data$time))^2)/n
    sigma2 = v0*tau2/chi
    sigma0 = solve(diag(8)*lambda)
    beta = rmvnorm(1, mean = mu_7, sigma = sigma2*sigma0)
    reg_curve[,i] = x %*% t(beta)
  }
  plot(data$time, Y, xlab="Time", ylab="Temperature in Celcius",
       main = TeX(paste0("$\\Omega_0 \ = \ $",lambda, "$\\I_3, \\v_0 \ = \ $",v0,
                         "$, \\sigma_0^2 \ = \ $",sigma20)))
  for (i in 1:size) {
    lines(data$time,reg_curve[,i],type = "l", col=rgb(0.9,0,0,0.2))
  }
}
prior_7(mu = mu_7)
