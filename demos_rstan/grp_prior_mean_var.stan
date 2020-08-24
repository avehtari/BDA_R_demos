
// Comparison of k groups with unequal variance and
// hierarchical priors for the mean and the variance
data {
  int<lower=0> N; // number of data points
  int<lower=0> K; // number of groups
  int<lower=1,upper=K> x[N]; // group indicator
  vector[N] y; //
}
parameters {
  real mu0;                 // prior mean
  real<lower=0> musigma0;   // prior std
  vector[K] mu;             // group means
  real lsigma0;             // prior mean
  real<lower=0> lsigma0s;   // prior std
  vector<lower=0>[K] sigma; // group stds
}
model {
  mu0 ~ normal(10, 10);        // weakly informative prior
  musigma0 ~ normal(0, 10);    // weakly informative prior
  mu ~ normal(mu0, musigma0);  // population prior with unknown parameters
  // lsigma0 is prior mean in log scale
  lsigma0 ~ normal(0, 1);      // weakly informative prior
  // log-normal prior sets normal prior on logarithm of the paremeter,
  // which is useful for positive parameters that shouldn't be very
  // close to 0. BDA3 Chapter 5 uses scaled inverse Chi^2 prior, but
  // as these don't need to be (semi-)conjugate thinking in terms of
  // log-normal can be easier.
  lsigma0s ~ lognormal(log(0.1), .5);   // weakly informative prior
  sigma ~ lognormal(lsigma0, lsigma0s); // population prior with unknown parameters
  y ~ normal(mu[x], sigma[x]); // observation model / likelihood
}
