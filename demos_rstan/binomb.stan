// Binomial model with a roughly uniform prior for
// the probability of success (theta)
data {
  int<lower=0> N;              // number of experiments
  int<lower=0> y;              // number of successes
}
parameters {
  // sampling is done for the parameters
  real alpha;                  // logit of probability of success in rang (-Inf,Inf)
}
transformed parameters {
  // transformed parameters are deterministic transformations of parameters (and data)
  real theta = inv_logit(alpha); // probability of success in range (0,1)
}
model {
  // model block creates the log density to be sampled
  alpha ~ normal(0, 1.5);        // roughly uniform prior for theta
  y ~ binomial_logit(N, alpha);  // model parameterized with logit of probability
  // the notation using ~ is syntactical sugar for
  //  target += normal_lpdf(alpha | 0, 1.5);       // lpdf for continuous theta
  //  target += binomial_logit_lpmf(y | N, alpha); // lpmf for discrete y
  // target is the log density to be sampled
}
