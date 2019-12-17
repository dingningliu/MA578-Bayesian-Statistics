data {
  int<lower=0> J;
  
  real y[J];
  
  real x[J];
}
parameters {
  vector[J] theta;
  real mu;
  real<lower=0> tau;
  real<lower=0> sigma;
  real a;
  real b;
}
model {
  for (j in 1:J)
     x[j] ~ normal(theta[j], sigma);
  for (j in 1:J)
     y[j] ~ normal(b*theta[j]+a, sqrt(2)*sigma);
  for (j in 1:J) theta[j] ~ normal(mu, tau);
  
}
