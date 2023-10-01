#' ---
#' title: "Bayesian data analysis demo 11.5"
#' author: "Aki Vehtari"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     theme: readable
#'     code_download: true
#' ---

#+ setup, message=FALSE, error=FALSE, warning=FALSE
library(posterior)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))

#' ## Use of posterior package for convergence diagnostics
#' 

#' In the Metropolis assignment it is likely that you have used a for
#' loop for the Metropolis algorithm. Here the Metropolis algorithm
#' has been replaced with independent draws from a normal
#' distribution.
#'
#' First illustration with a single chain
S <- 1000 # number of draws
D <- 2 # number of parameters
# 2-D matrix with dimensions ‘"iteration"’, and ‘"variable"’
thetas <- matrix(nrow=S, ncol=2)
colnames(thetas) <- c("theta1","theta2")
# Fill the matrix
for (s in 1:S) {
  thetas[s,] <- rnorm(n=2, mean=0, sd=1)
}

#' We can quickly get several useful summaries including
#' Rank-normalized-Rhat and ESS-bulk and ESS-tail
summarise_draws(thetas)

#' To compute the basic (split) Rhat and basic ESS as presented in the
#' lecture we can use additionl arguments
summarise_draws(thetas, Rhat=rhat_basic, ESS=ess_basic)

#' Compute mean and 5%- and 95%-quantiles and corresponding MCSE's
summarise_draws(thetas, mean, mcse_mean,
                ~quantile(.x, probs = c(0.05, 0.95)),
                ~mcse_quantile(.x, probs = c(0.05, 0.95)))

#' summarise_draws() function coerced the R matrix to a posterior
#' object behind the scenes, but we can also convert it explicitlu to
#' draws object that is better aware of number of iterations and
#' chains. For example, CmdStanR returns posterior draws objects.
#' posterior package provides several functions that make manipulating
#' posterior objects easier than plain R matrix and array types.
thetas <- as_draws_df(thetas)



#' Illustration with several chains
N <- 1000 # number of iterations per chain
M <- 4 # number of chains
D <- 2 # number of parameters
# 3-D array with dimensions ‘"iteration"’, ‘"chain"’, and ‘"variable"’
thetas <- array(dim=c(N,M,D))
# Fill the array
for (m in 1:M) { # loop over chains
  for (n in 1:N) { # loop over iterations
    thetas[n,m,] <- rnorm(n=2, mean=0, sd=1)
  }
}

#' We can quickly get several useful summaries including
#' Rank-normalized-Rhat and ESS-bulk and ESS-tail
summarise_draws(thetas)

#' To compute the basic (split) Rhat and basic ESS as presented in the
#' lecture we can use additionl arguments
summarise_draws(thetas, Rhat=rhat_basic, ESS=ess_basic)

#' Compute mean and 5%- and 95%-quantiles and corresponding MCSE's
summarise_draws(thetas, mean, mcse_mean,
                ~quantile(.x, probs = c(0.05, 0.95)),
                ~mcse_quantile(.x, probs = c(0.05, 0.95)))

#' summarise_draws() function coerced the R 3-D array to a posterior
#' object behind the scenes, but we can also convert it explicitlu to
#' draws object that is better aware of number of iterations and
#' chains. For example, CmdStanR returns posterior draws objects.
#' posterior package provides several functions that make manipulating
#' posterior objects easier than plain R matrix and array types.
thetas <- as_draws_df(thetas)
variables(thetas) <- c("theta1","theta2")

#' bayesplot package provides Some useful plots to analyse your MCMC
#' draws.  You don't need to include all of these in your assignment
#' or project reports, but they are shown as seeing the plo for your
#' own Metropolis algorithm helps to gain better understanding.
#'
#' Trace plots
#' 
mcmc_trace(thetas)
#'
#' Autocorrelation plots
#' 
mcmc_acf(thetas)
#'
#' Density plots with separate density for each chain
#' 
mcmc_dens_overlay(thetas)
#'
#' Density plots using all chains
#' 
mcmc_dens(thetas)
#'
#' Density plots with median and 50% central interval (see ?mcmc_areas
#' for more options)
#' 
mcmc_areas(thetas)
#' 
#' Scatter plots and histograms (works also for more than two variables)
#' 
mcmc_pairs(thetas)
