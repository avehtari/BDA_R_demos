// Gaussian linear model with standardized data
data {
  int<lower=0> N; // number of data points
  vector[N] x;    // covariate / predictor
  vector[N] y;    // target
  real xpred;     // new covariate value to make predictions
}
transformed data {
  // deterministic transformations of data
  vector[N] x_std = (x - mean(x)) / sd(x);
  vector[N] y_std = (y - mean(y)) / sd(y);
  real xpred_std  = (xpred - mean(x)) / sd(x);
}
parameters {
  real alpha;              // intercept
  real beta;               // slope
  real<lower=0> sigma_std; // standard deviation is constrained to be positive
}
transformed parameters {
  // deterministic transformation of parameters and data
  vector[N] mu_std = alpha + beta*x_std; // linear model
}
model {
  alpha ~ normal(0, 1); // weakly informative prior (given standardized data)
  beta ~ normal(0, 1);  // weakly informative prior (given standardized data)
  sigma_std ~ normal(0, 1); // weakly informative prior (given standardized data)
  y_std ~ normal(mu_std, sigma_std); // observation model / likelihood
}
generated quantities {
  // transform to the original data scale
  vector[N] mu = mu_std*sd(y) + mean(y);
  real<lower=0> sigma = sigma_std*sd(y);
  // sample from the predictive distribution
  real ypred = normal_rng((alpha + beta*xpred_std)*sd(y)+mean(y), sigma_std*sd(y));
  // compute log predictive densities to be used for LOO-CV
  // to make appropriate comparison to other models, this log density is computed
  // using the original data scale (y, mu, sigma)
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = normal_lpdf(y[i] | mu[i], sigma);
}
