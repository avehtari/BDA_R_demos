// Binomial model with beta(1,1) prior
data {
  int<lower=0> N;              // number of experiments
  int<lower=0> y;              // number of successes
}
parameters {
  real<lower=0,upper=1> theta; // probability of success in range (0,1)
}
model {
  // model block creates the log density to be sampled
  theta ~ beta(1, 1);          // prior
  y ~ binomial(N, theta);      // observation model / likelihood
  // the notation using ~ is syntactic sugar for
  //  target += beta_lpdf(theta | 1, 1);     // lpdf for continuous theta
  //  target += binomial_lpmf(y | N, theta); // lpmf for discrete y
  // target is the log density to be sampled
}
