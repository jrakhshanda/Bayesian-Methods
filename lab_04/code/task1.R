library(ggplot2)
library(rstan)
setwd("C:/Users/Rakshanda/Desktop/Bayesian Learning/LABS/lab_04")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#_____(a)_____#
library(latex2exp)
library(coda)

AR_1 = function(nDraws=200,mu=10,sigma2=2,phi){
  x = vector(length = nDraws)
  x[1] = mu
  for (i in 2:nDraws) {
    x[i] = mu + phi*(x[i-1]-mu) + rnorm(n = 1, mean = 0, sd = sqrt(sigma2))
  }
  ts.plot(x,ylab = 'x',
          main = TeX(paste0("$\\phi\ = \ $",phi)))
  return(x)
}

a = AR_1(phi = -0.95)
b = AR_1(phi = -0.5)
x = AR_1(phi = 0.3)
y = AR_1(phi = 0.95)
z = AR_1(phi = 1)
autocorr.plot(as.mcmc(a), main = 'Autocorrelation')
autocorr.plot(as.mcmc(b),main = 'Autocorrelation')
autocorr.plot(as.mcmc(x), main = 'Autocorrelation')
autocorr.plot(as.mcmc(y), main = 'Autocorrelation')
autocorr.plot(as.mcmc(z), main = 'Autocorrelation')
#______(b)_______#
# The input data is a vector 'y' of length 'N'
# The input data is a vector of length n
# MOdels accepts three unknown parameters

MCMC_posterior = '
data {
  int<lower=0> N;
  vector[N] X;
  
}

parameters {
  real mu;
  real<lower=0> sigma2;
  real phi;
  
}

model {
  for(i in 2:N) {
    X[i] ~  normal(mu + phi*(X[i-1] - mu), sigma2);
    
  }
}
'

x_t = list(N = length(x), X=x)
fit1 = stan(model_code = MCMC_posterior, data = x_t, 
            iter = 1000, warmup = 200, chains = 1)
library(coda)

# Print the fitted model
print(fit1,digits_summary=2)

# draws1 = extract(fit1) 
draws1 =As.mcmc.list(fit1)
plot(as.mcmc(draws1[,1:3]))

# joint posterior of mu and phi x_t
pairs(fit1,pars = c("mu","sigma2", "phi"))s

#____(ii)_____#
y_t = list(N = length(y), X=y)
fit2 = stan(model_code =  MCMC_posterior, data = y_t, 
            iter = 1000, warmup = 200, chains = 1)
print(fit2, digits_summary = 2)

#draws2 = extract(fit2) 
draws2 = As.mcmc.list(fit2)

# convergence of the sample x_t
plot(as.mcmc(draws2[,1:3]))

# joint posterior of mu and phi x_t

pairs(fit2,pars = c("mu", "phi"))

#_______(c)______#
data = read.table("campy.dat", header = T)
campy = as.vector(data$c)
N = length(campy)
postData = list(N = N, c = campy)

poisson_stan = '
data {
  int<lower=0> N;
  int  c[N];
}

parameters {
  real x[N];
  real mu;
  real phi;
  real<lower=0> sigma2;
  
}

model {
  
  for(i in 2:N){
  x[i] ~ normal(mu + phi*(x[i-1] - mu), sqrt(sigma2));
  }
  
  for(i in 1:N) {
    c[i] ~ poisson(exp(x[i]));
    
  }
}
'
fit = stan(model_code = poisson_stan, data = postData,
           iter = 1000, warmup = 500, chains = 1)
print(fit, digits_summary = 2)

sm = summary(fit)$summary
Mean = unlist(sapply(c(sm[1:N,1]), function(x) {exp(x)}))
lower = unlist(sapply(c(sm[1:N,4]), function(x) {exp(x)}))
upper = unlist(sapply(c(sm[1:N,8]), function(x) {exp(x)}))


df = cbind(index = 1:N,campy, Mean, lower, upper)
ggplot(as.data.frame(df))+ theme_classic()+
  geom_point(mapping = aes(x = index, y = campy))+
  geom_line(mapping = aes(x = index, y = Mean, col = 'Mean'), size = 1)+
  geom_ribbon(mapping = aes(x= index, ymin= lower, ymax =upper),
              alpha = 0.4)
#________(d)_________#
poissonPosterior_prior = '

data {
  int<lower=0> N;
  int  c[N];
  real<lower=0> v0;
  real<lower=0> sigma20;
}

parameters {
  real x[N];
  real mu;
  real phi;
  real<lower=0> sigma2;
  
}

model {

  sigma2 ~ scaled_inv_chi_square(v0, sqrt(sigma20));
  
  for(i in 2:N){
  x[i] ~ normal(mu + phi*(x[i-1] - mu), sqrt(sigma2));
  }
  
  for(i in 1:N) {
    c[i] ~ poisson(exp(x[i]));
    
  }
}
'
priorData = list(N = N, c = campy, v0 = 50, sigma20 = 0.01)

fit1 = stan(model_code = poissonPosterior_prior,data = priorData, 
            iter = 1000, warmup = 500)

print(fit1)
sm1 = summary(fit1)$summary
Mean1 = unlist(sapply(c(sm1[1:N,1]), function(x) {exp(x)}))
lower1 = unlist(sapply(c(sm1[1:N,4]), function(x) {exp(x)}))
upper1 = unlist(sapply(c(sm1[1:N,8]), function(x) {exp(x)}))


df = cbind(index = 1:N,campy, Mean1, lower1, upper1)
ggplot(as.data.frame(df))+ theme_classic()+
  geom_point(mapping = aes(x = index, y = campy))+
  geom_line(mapping = aes(x = index, y = Mean1, col = 'Mean'), size = 1)+
  geom_ribbon(mapping = aes(x= index, ymin= lower1, ymax =upper1),
              alpha = 0.4)



