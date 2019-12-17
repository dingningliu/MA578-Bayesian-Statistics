data {
  int<lower=0> J;
  
  real y[J];
  
  real x[J];
  
  real z[J];
}
parameters {
  vector[J] u;
  vector[J] v;
  real mu1;
  real mu2;
  real<lower=0> tau1;
  real<lower=0> tau2;
  real<lower=0> sigma;
  real a;
  real b;
  real c;
}
model {
  for (j in 1:J)
     x[j] ~ normal(u[j], sigma);
  for (j in 1:J)
     z[j] ~ normal(v[j], sigma);
  for (j in 1:J)
     y[j] ~ normal(b*u[j]+c*v[j]+a, sigma);
  for (j in 1:J) u[j] ~ normal(mu1, tau1);
  for (j in 1:J) v[j] ~ normal(mu2, tau2);
  
}
