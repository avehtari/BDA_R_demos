//  Comparison of two groups with Binomial
data {
  int<lower=0> N1;              // number of experiments in group 1
  int<lower=0> y1;              // number of deaths in group 1
  int<lower=0> N2;              // number of experiments in group 2
  int<lower=0> y2;              // number of deaths in group 2
}
parameters {
  real<lower=0,upper=1> theta1; // probability of death in group 1
  real<lower=0,upper=1> theta2; // probability of death in group 2
}
model {
  // model block creates the log density to be sampled
  theta1 ~ beta(1, 1);          // prior
  theta2 ~ beta(1, 1);          // prior
  y1 ~ binomial(N1, theta1);    // observation model / likelihood
  y2 ~ binomial(N2, theta2);    // observation model / likelihood
  // the notation using ~ is syntactical sugar for
  //  target += beta_lpdf(theta1 | 1, 1);       // lpdf for continuous theta1
  //  target += beta_lpdf(theta2 | 1, 1);       // lpdf for continuous theta2
  //  target += binomial_lpmf(y1 | N1, theta1); // lpmf for discrete y1
  //  target += binomial_lpmf(y2 | N2, theta2); // lpmf for discrete y2
  // target is the log density to be sampled
}
generated quantities {
  // generated quantities are computed after sampling
  real oddsratio = (theta2/(1-theta2))/(theta1/(1-theta1));
}
