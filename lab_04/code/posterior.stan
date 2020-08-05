
// The input data is a vector 'X' of length 'N'

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
    X[i] ~  normal(mu + phi*(X[i-1] - mu), sqrt(sigma2));
    
  }
}