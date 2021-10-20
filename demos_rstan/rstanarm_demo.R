#' ---
#' title: "Bayesian data analysis - RStanARM demos"
#' author: "Aki Vehtari, Markus Paasiniemi"
#' date: "First version 2017-07-17. Last modified `r format(Sys.Date())`."
#' output:
#'   html_document:
#'     fig_caption: yes
#'     toc: TRUE
#'     toc_depth: 2
#'     number_sections: TRUE
#'     toc_float:
#'       smooth_scroll: FALSE
#'     theme: readable
#'     code_download: true
#' ---

#' # Setup  {.unnumbered}

#+ setup, include=FALSE
knitr::opts_chunk$set(cache=FALSE, message=FALSE, error=FALSE, warning=TRUE, comment=NA, out.width='95%')


#' **Load packages**
library(tidyr)
library(dplyr)
library(rstan) 
library(rstanarm)
options(mc.cores = 1)
library(loo)
library(shinystan)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(ggdist)
library(gridExtra)
library(rprojroot)
root<-has_file(".BDA_R_demos_root")$make_fix_file()
SEED <- 48927 # set random seed for reproducability


#' # Introduction
#' 
#' This notebook contains several examples of how to use [Stan](https://mc-stan.org) in R with __rstanarm__. This notebook assumes basic knowledge of Bayesian inference and MCMC. The examples are related to [Bayesian data analysis course](https://avehtari.github.io/BDA_course_Aalto/).
#' 
#' Note that you can easily analyse Stan fit objects returned by `stan_glm()` with a ShinyStan package by calling `launch_shinystan(fit)`.
#' 
#' The models are not exactly equal to the models at rstan_demo.Rmd, but rather serve as examples of how to implement similar models with __rstanarm__.
#' 
#' # Bernoulli model
#' 
#' Toy data with sequence of failures (0) and successes (1). We would like to learn about the unknown probability of success.
data_bern <- data.frame(y = c(1, 1, 1, 0, 1, 1, 1, 0, 1, 0))


#' Uniform prior (beta(1,1)) is achieved by setting the prior to NULL,
#' which is not recommended in general. y ~ 1 means y depends only on
#' the intercept term
fit_bern <- stan_glm(y ~ 1, family = binomial(), data = data_bern,
                     prior_intercet = NULL, seed = SEED, refresh = 0)


#' You can use ShinyStan examine and diagnose the fitted model is to call shinystan in R terminal as `launch_shinystan(fit_bern)`

#' Monitor provides summary statistics and diagnostics
monitor(fit_bern$stanfit)


#' To see the parameter values on the ouput space, do the inverse
#' logistic transformation (plogis in R) on the intercept
draws <- as.data.frame(fit_bern)
mean(draws$`(Intercept)`)

#' Probability of success
draws$theta <- plogis(draws$`(Intercept)`)
mean(draws$theta)


#' Histogram of theta
mcmc_hist(draws, pars='theta') + xlab('theta')


#' We next compare the result to using the default prior which is normal(0, 2.5) on logit probability. Visualize the prior by drawing samples from it
prior_mean <- 0
prior_sd <- 2.5
prior_intercept <- normal(location = prior_mean, scale = prior_sd)
prior_samples <- data.frame(
                 theta = plogis(rnorm(20000, prior_mean, prior_sd)))
mcmc_hist(prior_samples)


fit_bern <- stan_glm(y ~ 1, family = binomial(), data = data_bern,
                     seed = SEED, refresh = 0)


monitor(fit_bern$stanfit)


#' To see the parameter values on the ouput space, do the inverse
#' logistic transformation (plogis in R) on the intercept
draws <- as.data.frame(fit_bern)
mean(draws$`(Intercept)`)

#' Probability of success
draws$theta <- plogis(draws$`(Intercept)`)
mean(draws$theta)


#' Histogram of theta
mcmc_hist(draws, pars='theta') + xlab('theta')

#' As the number of observations is small, there is small change in the posterior mean when the prior is changed. You can experiment with different priors and varying the number of observations.
#' 
#' 
#' # Binomial model
#' 
#' Instead of sequence of 0's and 1's, we can summarize the data with the number of experiments and the number successes. Binomial model with a approximately uniform prior for the probability of success. The prior is specified in the 'latent space'. The actual probability of success, theta = plogis(alpha), where plogis is the inverse of the logistic function.
#' 
#' Visualize the prior by drawing samples from it
prior_mean <- 0
prior_sd <- 1.5
prior_intercept <- normal(location = prior_mean, scale = prior_sd)
prior_samples <- data.frame(
                 theta = plogis(rnorm(20000, prior_mean, prior_sd)))
mcmc_hist(prior_samples)


#' Binomial model (we are not able to replicate the Binomial example in rstan_demo exactly, as `stan_glm` does not accept just one observation, so the Bernoulli is needed for the same model, and Binomial will be demonstrated first with other data).
data_bin <- data.frame(N = c(5,5), y = c(4,3))
fit_bin <- stan_glm(y/N ~ 1, family = binomial(), data = data_bin,
                     prior_intercept = prior_intercept, weights = N,
		     seed = SEED, refresh = 0)


monitor(fit_bin$stanfit)


draws <- as.data.frame(fit_bin)
mean(draws$`(Intercept)`)

#' Probability of success
draws$theta <- plogis(draws$`(Intercept)`)
mean(draws$theta)


#' Histogram of theta
mcmc_hist(draws, pars='theta') + xlab('theta')


#' Re-run the model with a new data dataset.
data_bin <- data.frame(N = c(5,5), y = c(4,5))
fit_bin <- update(fit_bin, data = data_bin)

monitor(fit_bin$stanfit)


#' Probability of success
draws <- as.data.frame(fit_bern)
draws$theta <- plogis(draws$`(Intercept)`)
mean(draws$theta)


#' Histogram of theta
mcmc_hist(draws, pars='theta') + xlab('theta')


#' # Comparison of two groups with Binomial 
#' 
#' An experiment was performed to estimate the effect of beta-blockers on mortality of cardiac patients. A group of patients were randomly assigned to treatment and control groups:
#' 
#' - out of 674 patients receiving the control, 39 died
#' - out of 680 receiving the treatment, 22 died
#' 
#' Data, where grp2 is a dummy variable that captures the differece of
#' the intercepts in the first and the second group.
data_bin2 <- data.frame(N = c(674, 680), y = c(39,22), grp2 = c(0,1))


#' To analyse whether the treatment is useful, we can use Binomial model for both groups and compute odds-ratio.

fit_bin2 <- stan_glm(y/N ~ grp2, family = binomial(), data = data_bin2,
                     weights = N, seed = SEED, refresh = 0)

monitor(fit_bin2$stanfit)


#' Plot odds ratio
draws_bin2 <- as.data.frame(fit_bin2) %>%
  mutate(theta1 = plogis(`(Intercept)`),
         theta2 = plogis(`(Intercept)` + grp2),
         oddsratio = (theta2/(1-theta2))/(theta1/(1-theta1)))
mcmc_hist(draws_bin2, pars='oddsratio')


#' # Linear Gaussian model
#' 
#' The following file has Kilpisjärvi summer month temperatures 1952-2013:
data_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
data_lin <-data.frame(year = data_kilpis$year,
                   temp = data_kilpis[,5])


#' Plot the data
ggplot() +
  geom_point(aes(year, temp), data = data.frame(data_lin), size = 1) +
  labs(y = 'Summer temp. @Kilpisjärvi', x= "Year") +
  guides(linetype = "none")


#' To analyse has there been change in the average summer month temperature we use a linear model with Gaussian model for the unexplained variation. rstanarm uses by default scaled priors.
#' 
#' y ~ x means y depends on the intercept and x
fit_lin <- stan_glm(temp ~ year, data = data_lin, family = gaussian(),
                    seed = SEED, refresh = 0)


#' The default priors for the linear model are
prior_summary(fit_lin)


#' 
#' You can use ShinyStan (`launch_shinystan(fit_lin)`) to look at the divergences, treedepth exceedences, n_eff, Rhats, and joint posterior of alpha and beta. In the corresponding rstan_demo notebook we observed some treedepth exceedences leading to slightly less efficient sampling, but rstanarm has slightly different model and performs better.
#' 
#' Instead of interactive ShinyStan, we can also check the diagnostics as follows
monitor(fit_lin$stanfit)
check_hmc_diagnostics(fit_lin$stanfit)


#' Plot data and the fit
draws_lin <- as.data.frame(fit_lin)
mean(draws_lin$year>0) # probability that beta > 0

mu_draws <- tcrossprod(cbind(1, data_lin$year),
                         cbind(draws_lin$`(Intercept)`,draws_lin$year))
mu <- apply(mu_draws, 1, quantile, c(0.05, 0.5, 0.95)) %>%
  t() %>% data.frame(x = data_lin$year, .) %>% gather(pct, y, -x)
pfit <- ggplot() +
  geom_point(aes(year, temp), data = data.frame(data_lin), size = 1) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(x = '', y = 'Summer temp. @Kilpisjärvi') +
  guides(linetype = "none")
phist <- mcmc_hist(draws_lin) + ggtitle('parameters')
grid.arrange(pfit, phist)


#' Prediction for year 2016
predict(fit_lin, newdata = data.frame(year = 2016), se.fit = TRUE)
# or sample from the posterior predictive distribution and
# plot the histogram
ypred <- posterior_predict(fit_lin, newdata = data.frame(year = 2016))
mcmc_hist(ypred) + xlab('avg-temperature prediction for the summer 2016')


#' # Linear Student's t model with brms
#' 
#' The temperatures used in the above analyses are averages over three months, which makes it more likely that they are normally distributed, but there can be extreme events in the feather and we can check whether more robust Student's t observation model woul give different results.
#' 
#' Currently, rstanarm does not yet support Student's t likelihood. Below we use brms package, which supports similar model formulas as rstanarm with more options, but doesn't have pre-compiled models (be aware also that the default priors are not necessary the same).

#+  results='hide'
library(brms)
fit_lin_t <- brm(temp ~ year, data = data_lin, family = student(), seed = SEED,
                 refresh = 1000)


summary(fit_lin_t)


#' brms package generates Stan code which we can extract as follows. By saving this code to a file you can extend the model, beyond the models supported by brms.
stancode(fit_lin_t)


#' # Pareto-smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO)
#' 
#' We can use leave-one-out cross-validation to compare the expected predictive performance.
#' 
#' Let's use LOO to compare whether Student's t model has better predictive performance.
loo1 <- loo(fit_lin)
loo2 <- loo(fit_lin_t)
loo_compare(loo1, loo2)

#' There is no practical difference between Gaussian and Student's t models.
#' 
#' # Comparison of k groups with hierarchical models
#' 
#' Let's compare the temperatures in three summer months.
data_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
data_grp <- data.frame(month = rep(6:8, nrow(data_kilpis)),
              temp = c(t(data_kilpis[,2:4])))


#' # Common variance (ANOVA) model
#' 
#' Weakly informative prior for the common mean
prior_intercept <- normal(10, 10)

#' To use no (= uniform) prior, prior_intercept could be set to NULL
#' 
#' y ~ 1 + (1 | x) means y depends on common intercept and group speficific intercepts (grouping determined by x)
fit_grp <- stan_lmer(temp ~ 1 + (1 | month), data = data_grp,
                     prior_intercept = prior_intercept, refresh = 0)
# launch_shinystan(fit_grp)


monitor(fit_grp$stanfit)

#' Average temperature over all months. monthly deviations from the
#' mean, residual sigma and hierarchical prior sigma
mcmc_hist(as.data.frame(fit_grp))

#' A density estimates of the posterior for each month
temps <- (as.matrix(fit_grp)[,1] + as.matrix(fit_grp)[, 2:4]) %>%
  as.data.frame() %>% setNames(c('June','July','August')) %>% gather(month, temp)
ggplot(temps, aes(y=month, x=temp)) +
  stat_slab() + labs(y='Month', x='Temperature')

#' Probabilities that June is hotter than July, June is hotter than August
#' and July is hotter than August:
combn(unique(temps$month), 2, function(months, data) {
  mean(subset(data, month == months[1])$temp > subset(data, month == months[2])$temp)
}, data = temps) %>% setNames(c('TJune>TJuly', 'TJune>TAugust', 'TJuly>TAugust'))


#' <br />
#' 
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2017-2019, Aki Vehtari, 2017 Markus Paasiniemi, licensed under BSD-3.
#' * Text &copy; 2017-2019, Aki Vehtari, licensed under CC-BY-NC 4.0.
