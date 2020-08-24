// Bernoulli model
data {
  int<lower=0> N;              // number of observations
  int<lower=0,upper=1> y[N];   // vector of binary observations
}
parameters {
  real<lower=0,upper=1> theta; // probability of success
}
model {
  // model block creates the log density to be sampled
  theta ~ beta(1, 1);          // prior 
  y ~ bernoulli(theta);        // observation model / likelihood
  // the notation using ~ is syntactical sugar for
  //  target += beta_lpdf(theta | 1, 1);   // lpdf for continuous theta
  //  target += bernoulli_lpmf(y | theta); // lpmf for discrete y
  // target is the log density to be sampled
}
