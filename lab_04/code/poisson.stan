data {
  int<lower=0> N;
  int  y[N];
  
}

parameters {
  real x[N];
  real mu;
  real phi;
  real sigma2;
  
}

model {
  x[1] ~ normal(mu, sigma2);
  y[1] ~ poisson(exp(x[1]));
  for(i in 2:N) {
    x[i] ~ normal(mu + phi*(x[i-1] - mu), sigma2);
    y[i] ~ poisson(exp(x[i]));
    
  }
}