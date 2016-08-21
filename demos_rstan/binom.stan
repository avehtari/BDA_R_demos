// Binomial model with a roughly uniform prior for
// the probability of success (theta)
data {
  int<lower=0> N;
  int<lower=0> y;
}
parameters {
  real alpha;
}
transformed parameters {
  real theta;
  theta = inv_logit(alpha);
}
model {
  // roughly auniform prior for the number of successes
  // see r-code for a visualization of the prior
  alpha ~ normal(0,1.5);
  y ~ binomial_logit(N,alpha);
}
