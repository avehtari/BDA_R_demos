// Gaussian linear model with adjustable priors
data {
  int<lower=0> N; // number of data points
  vector[N] x; // covariate / predictor
  vector[N] y; // target
  real xpred; // new covariate value to make predictions
  real pmualpha; // prior mean for alpha
  real psalpha; // prior std for alpha
  real pmubeta; // prior mean for beta
  real psbeta; // prior std for beta
  real pssigma; // prior std for half-normal prior for sigma
}
parameters {
  real alpha; // intercept
  real beta; // slope
  real<lower=0> sigma; // standard deviation is constrained to be positive
}
transformed parameters {
  // deterministic transformation of parameters and data
  vector[N] mu = alpha + beta * x; // linear model
}
model {
  alpha ~ normal(pmualpha, psalpha); // prior
  beta ~ normal(pmubeta, psbeta); // prior
  sigma ~ normal(0, pssigma); // as sigma is constrained to be positive,
  // this is same as half-normal prior
  y ~ normal(mu, sigma); // observation model / likelihood
  // the notation using ~ is syntactic sugar for
  //  target += normal_lpdf(alpha | pmualpha, psalpha);
  //  target += normal_lpdf(beta | pmubeta, psbeta);
  //  target += normal_lpdf(y | mu, sigma);
  // target is the log density to be sampled
}
generated quantities {
  // sample from the predictive distribution
  real ypred = normal_rng(alpha + beta * xpred, sigma);
  // compute log predictive densities to be used for LOO-CV
  vector[N] log_lik;
  for (i in 1 : N) {
    log_lik[i] = normal_lpdf(y[i] | mu[i], sigma);
  }
}
