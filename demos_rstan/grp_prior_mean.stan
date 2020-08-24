// Comparison of k groups with common variance and
// hierarchical prior for the mean
data {
  int<lower=0> N;            // number of observations
  int<lower=0> K;            // number of groups
  int<lower=1,upper=K> x[N]; // discrete group indicators
  vector[N] y;               // real valued observations
}
parameters {
  real mu0;             // prior mean
  real<lower=0> sigma0; // prior std constrained to be positive
  vector[K] mu; // group means
  real<lower=0> sigma;  // common std constrained to be positive
}
model {
  mu0 ~ normal(10, 10);     // weakly informative prior
  sigma0 ~ normal(0, 10);   // weakly informative prior
  mu ~ normal(mu0, sigma0); // population prior with unknown parameters
  // log-normal prior sets normal prior on logarithm of the paremeter,
  // which is useful for positive parameters that shouldn't be very
  // close to 0. BDA3 Chapter 5 uses scaled inverse Chi^2 prior, but
  // as these don't need to be (semi-)conjugate, thinking in terms of
  // log-normal can be easier.
  sigma ~ lognormal(0, .5); // weakly informative prior
  y ~ normal(mu[x], sigma); // observation model / likelihood
}
