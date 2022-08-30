// Comparison of k groups with common variance (ANOVA)
data {
  int<lower=0> N; // number of observations
  int<lower=0> K; // number of groups
  array[N] int<lower=1, upper=K> x; // discrete group indicators
  vector[N] y; // real valued observations
}
parameters {
  vector[K] mu;        // group means
  real<lower=0> sigma; // common standard deviation constrained to be positive
}
model {
  mu ~ normal(0, 100); // weakly informative prior
  sigma ~ normal(0, 1); // weakly informative prior
  y ~ normal(mu[x], sigma); // observation model / likelihood
}
