#' ---
#' title: "Bayesian data analysis - CmdStanR demos"
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
library(cmdstanr)
# If running in Aalto JupyterHub the path should automatically be set to
# '/coursedata/cmdstan'
options(mc.cores = 1)
library(posterior)
options(posterior.num_args=list(sigfig=2)) # by default summaries with 2 significant digits
library(loo)
library(tidyr)
library(dplyr)
options(pillar.neg=FALSE)
library(ggplot2)
library(gridExtra)
library(bayesplot)
library(ggdist)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(rprojroot)
root<-has_file(".BDA_R_demos_root")$make_fix_file()
SEED <- 48927 # set random seed for reproducability


#' # Introduction
#' 
#' This notebook contains several examples of how to use [Stan](https://mc-stan.org) in R with __cmdstanr__. This notebook assumes basic knowledge of Bayesian inference and MCMC. The Stan models are stored in separate .stan-files. The examples are related to [Bayesian data analysis course](https://avehtari.github.io/BDA_course_Aalto/).
#' 
#' # Bernoulli model
#' 
#' Toy data with sequence of failures (0) and successes (1). We would like to learn about the unknown probability of success.
data_bern <- list(N = 10, y = c(1, 1, 1, 0, 1, 1, 1, 0, 1, 0))


#' Bernoulli model with a proper Beta(1,1) (uniform) prior
code_bern <- root("demos_rstan", "bern.stan")
writeLines(readLines(code_bern))


#' Sample form the posterior and show the summary
#+ results='hide'
mod_bern <- cmdstan_model(stan_file = code_bern)
fit_bern <- mod_bern$sample(data = data_bern, seed = SEED, refresh=1000)

fit_bern$summary()

#' Plot a histogram of the posterior draws with bayesplot (uses ggplot)
draws <- fit_bern$draws(format = "df")
mcmc_hist(draws, pars='theta') + xlim(c(0,1))

#' Plot a dots plot of the posterior draws with ggplot + ggdist
#+ warning=FALSE
draws |>
  ggplot(aes(x=theta)) + 
  stat_dotsinterval() + 
  xlim(c(0,1))

#' # Binomial model
#' 
#' Instead of sequence of 0's and 1's, we can summarize the data with the number of experiments and the number successes:
data_bin <- list(N = 10, y = 7)


#' And then we use Binomial model with Beta(1,1) prior for the probability of success.
code_binom <- root("demos_rstan","binom.stan")
writeLines(readLines(code_binom))


#' Sample from the posterior and plot the posterior. The histogram should look similar as in the Bernoulli case.
#' 
#+ results='hide'
mod_bin <- cmdstan_model(stan_file = code_binom)
fit_bin <- mod_bin$sample(data = data_bin, seed = SEED, refresh=1000)

fit_bin$summary()

draws <- fit_bin$draws(format = "df")
mcmc_hist(draws, pars = 'theta') + xlim(c(0,1))


#' Re-run the model with a new data. The compiled Stan program is re-used making the re-use faster.
#' 
#+ results='hide'
data_bin <- list(N = 100, y = 70)
fit_bin <- mod_bin$sample(data = data_bin, seed = SEED, refresh=1000)

fit_bin$summary()

draws <- fit_bin$draws(format = "df")
mcmc_hist(draws, pars = 'theta', binwidth = 0.01) + xlim(c(0,1))


#' ## Explicit transformation of variables
#' 
#' In the above examples the probability of success $\theta$ was declared as
#' 
#' `real<lower=0,upper=1> theta;`
#' 
#' Stan makes automatic transformation of the variable to the unconstrained space using logit transofrmation for interval constrained and log transformation for half constraints.
#' 
#' The following example shows how we can also make an explicit transformation and use binomial_logit function which takes the unconstrained parameter as an argument and uses logit transformation internally. This form can be useful for better numerical stability.

code_binomb <- root("demos_rstan", "binomb.stan")
writeLines(readLines(code_binomb))

#' Here we have used Gaussian prior in the unconstrained space, which produces close to uniform prior for theta.
#' 
#' Sample from the posterior and plot the posterior. The histogram should look similar as with the previous models.
#' 
#+ results='hide'
data_bin <- list(N = 100, y = 70)
mod_binb <- cmdstan_model(stan_file = code_binomb)
fit_binb <- mod_bin$sample(data = data_bin, seed = SEED, refresh=1000)

fit_binb$summary()

draws <- fit_binb$draws(format = "df")
mcmc_hist(draws, pars = 'theta', binwidth = 0.01) + xlim(c(0,1))


#' 
#' # Comparison of two groups with Binomial
#' 
#' An experiment was performed to estimate the effect of beta-blockers on mortality of cardiac patients. A group of patients were randomly assigned to treatment and control groups:
#' 
#' - out of 674 patients receiving the control, 39 died
#' - out of 680 receiving the treatment, 22 died
#' 
#' Data:

data_bin2 <- list(N1 = 674, y1 = 39, N2 = 680, y2 = 22)


#' To analyse whether the treatment is useful, we can use Binomial model for both groups and compute odds-ratio:
code_binom2 <- root("demos_rstan", "binom2.stan")
writeLines(readLines(code_binom2))


#' Sample from the posterior and plot the posterior
#+ results='hide'
mod_bin2 <- cmdstan_model(stan_file = code_binom2)
fit_bin2 <- mod_bin2$sample(data = data_bin2, seed = SEED, refresh=1000)

fit_bin2$summary()

#' Histogram
#+  warning=FALSE
draws <- fit_bin2$draws(format = "df")
mcmc_hist(draws, pars = 'oddsratio') +
  geom_vline(xintercept = 1) +
  scale_x_continuous(breaks = c(seq(0.25,1.5,by=0.25)))

#' Dots plot with median and 66% and 95% intervals
#+ warning=FALSE
draws |>
  ggplot(aes(x=oddsratio)) + 
    geom_dotsinterval() + 
    geom_vline(xintercept = 1) +
    scale_x_continuous(breaks = c(seq(0.25,1.5,by=0.25)))+
    labs(x='Odds ratio', y='')

#' Probability (and corresponding MCSE) that oddsratio<1
draws |>
  mutate_variables(p_oddsratio_lt_1 = as.numeric(oddsratio<1)) |>
  subset_draws("p_oddsratio_lt_1") |>
  summarise_draws(prob=mean, MCSE=mcse_mean)

#' # Linear Gaussian model
#' 
#' The following file has Kilpisjärvi summer month temperatures 1952-2022
#' (data by Finnish Meteorological Institute, CC-BY 4.0). 
data_kilpis <- read.delim(root("demos_rstan","kilpisjarvi-summer-temp-2022.csv"), sep = ";")
data_lin <-list(N = nrow(data_kilpis),
             x = data_kilpis$year,
             xpred = 2016,
             y = data_kilpis[,5])


#' Plot the data
ggplot() +
  geom_point(aes(x, y), data = data.frame(data_lin), size = 1) +
  labs(y = 'Summer temp. @Kilpisjärvi', x= "Year") +
  guides(linetype = "none")


#' To analyse whether the average summer month temperature is rising, we use a linear model with Gaussian model for the unexplained variation. 
#' 
#' ## Gaussian linear model with adjustable priors
#' 
#' The folloing Stan code allows also setting hyperparameter values as data allowing easier way to use different priors in different analyses:
code_lin <- root("demos_rstan", "lin.stan")
writeLines(readLines(code_lin))

#' Create another list with data and priors
data_lin_priors <- c(list(
    pmualpha = mean(unlist(data_kilpis[,5])), # centered
    psalpha = 100, # weakly informative
    pmubeta = 0, # a priori incr. and decr. as likely
    psbeta = (.1--.1)/6, # avg temp prob does does not incr. more than a degree per 10 years
    pssigma = 1), # total variation in summer average temperatures is less +-3 degrees
  data_lin)


#' Run Stan
#+ results='hide'
mod_lin <- cmdstan_model(stan_file = code_lin)
fit_lin <- mod_lin$sample(data = data_lin_priors, seed = SEED, refresh=1000)

#' Stan gives a warning: There were X transitions after warmup that exceeded the maximum treedepth. 
#' 
#' We can check the generic convergence diagnostics as follows
fit_lin$summary()

#' We can check the HMC specific diagnostics as follows
#+  message=TRUE
fit_lin$diagnostic_summary(diagnostics = c("divergences", "treedepth"))

#' The high number of max treedepth exceedence is due to high posterior dependency,
#' which in this case reduces efficiency, but doesn't invalidate the results
draws_lin <- fit_lin$draws(format = "df")
draws_lin |>
  mcmc_scatter(pars=c("alpha","beta"))

#' 
#' Compute the probability that the summer temperature is increasing.
#+ warning=FALSE
mean(draws_lin[,"beta"]>0) # probability that beta > 0


#' Plot the data, the model fit and prediction for year 2016.
#+ warning=FALSE
mu <- draws_lin |>
  as_draws_df() |>
  as_tibble() |>
  select(starts_with("mu")) |>
  apply(2, quantile, c(0.05, 0.5, 0.95)) |>
  t() %>%
  data.frame(x = data_lin$x, .)  |> 
  gather(pct, y, -x)

pfit <- ggplot() +
  geom_point(aes(x, y), data = data.frame(data_lin), size = 1) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(y = 'Summer temp. @Kilpisjärvi', x= "Year") +
  guides(linetype = "none")
phist <- mcmc_hist(draws_lin, pars = c('beta','sigma','ypred'))
grid.arrange(pfit, phist, nrow = 2)


#' ## Gaussian linear model with standardized data
#' 
#' In the above we used the unnormalized data and as x values are far away from zero, this will cause very strong posterior dependency between alpha and beta (did you use ShinyStan for the above model?). The strong posterior dependency can be removed by normalizing the data to have zero mean. The following Stan code makes it in Stan. In generated quantities we do correspnding transformation back to the original scale.

code_lin_std <- root("demos_rstan", "lin_std.stan")
writeLines(readLines(code_lin_std))


#+ results='hide'
mod_lin_std <- cmdstan_model(stan_file = code_lin_std)
fit_lin_std <- mod_lin_std$sample(data = data_lin, seed = SEED, refresh=1000)

#' Now there were no warnings. We can check diagnostics with the following commands.

#+  message=TRUE
fit_lin_std$summary()
fit_lin_std$diagnostic_summary(diagnostics = c("divergences", "treedepth"))


#' We see that there are no warnings by diagnostics and ESS's are higher than with the previous case with non-standardized data. The posterior has no high dependency.
#'
draws_lin_std <- fit_lin_std$draws(format = "df")
draws_lin_std |>
  mcmc_scatter(pars=c("alpha","beta"))

#' 
#' Next we check that we get similar probability for beta>0.
#+ warning=FALSE
draws_lin_std <- fit_lin_std$draws(format = "df")
mean(draws_lin_std[,"beta"]>0) # probability that beta > 0


#' # Linear Student's $t$ model.
#' 
#' The temperatures used in the above analyses are averages over three months, which makes it more likely that they are normally distributed, but there can be extreme events in the feather and we can check whether more robust Student's $t$ observation model woul give different results.

code_lin_std_t <- root("demos_rstan", "lin_std_t.stan")
writeLines(readLines(code_lin_std_t))


#+ results='hide'
mod_lin_std_t <- cmdstan_model(stan_file = code_lin_std_t)
fit_lin_std_t <- mod_lin_std_t$sample(data = data_lin, seed = SEED, refresh=1000)

#' We get some warnings, but these specific warnings are not critical if counts are small as here.
#' 
#' Let's examine further diagnostics.
#+  message=TRUE
fit_lin_std_t$summary()
fit_lin_std_t$diagnostic_summary(diagnostics = c("divergences", "treedepth"))


#' We get similar diagnostics as for the linear Gaussian model with non-standardised data.
#' 
#' Compute the probability that the summer temperature is increasing.
#+ warning=FALSE
draws_lin_std_t <- fit_lin_std_t$draws(format = "df")
mean(draws_lin_std_t[,"beta"]>0) # probability that beta > 0

#' We get similar probability as with Gaussian obervation model.
#' 
#' 
#' Plot data and the model fit
#+ warning=FALSE
mu <- draws_lin_std_t |>
  as_tibble() |>
  select(starts_with("mu[")) |>
  apply(2, quantile, c(0.05, 0.5, 0.95)) |>
  t() |> 
  data.frame(x = data_lin$x, .)  |> 
  gather(pct, y, -x)

pfit <- ggplot() +
  geom_point(aes(x, y), data = data.frame(data_lin), size = 1) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(y = 'Summer temp. @Kilpisjärvi', x= "Year") +
  guides(linetype = "none")
phist <- mcmc_hist(draws_lin_std_t, pars = c('beta','sigma','ypred'))
grid.arrange(pfit, phist, nrow = 2)

#' We see also that the marginal posterior of nu is wide with lot of mass for values producing distrbution really close to Gaussian.
#' 
#' # Pareto-smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO)
#' 
#' We can use leave-one-out cross-validation to compare the expected predictive performance. For the following lines to work, the log-likelihood needs to be evaluated in the stan code. For an example, see lin.stan and [Computing approximate leave-one-out cross-validation usig PSIS-LOO](http://mc-stan.org/loo/articles/loo2-with-rstan.html).
# At this moment there is not yet method for cmdstanr fit, so we use this helper
loocmd <- function(fit, ...) {
  loo(fit$draws("log_lik"), r_eff=relative_eff(fit$draws("log_lik")), ...)
}
loo_lin_std <- loocmd(fit_lin_std)
loo_lin_std_t <- loocmd(fit_lin_std_t)
loo_compare(loo_lin_std, loo_lin_std_t)

#' There is no practical difference between Gaussian and Student's $t$ observation model for this data.
#' 
#' 
#' # Comparison of $k$ groups with hierarchical models
#' 
#' Let's compare the temperatures in three summer months.
data_grp <-list(N = 3*nrow(data_kilpis),
             K = 3,
             x = rep(1:3, nrow(data_kilpis)),
             y = c(t(data_kilpis[,2:4])))


#' ## Common variance (ANOVA) model
code_grp_aov <- root("demos_rstan", "grp_aov.stan")
writeLines(readLines(code_grp_aov))


#' Fit the model
#+ results='hide'
mod_grp <- cmdstan_model(stan_file = code_grp_aov)
fit_grp <- mod_grp$sample(data = data_grp, seed = SEED, refresh=1000)

fit_grp$summary()
fit_grp$diagnostic_summary(diagnostics = c("divergences", "treedepth"))

#' ## Common variance and hierarchical prior for mean.
#' 
#' Results do not differ much from the previous, because there is only
#' few groups and quite much data per group, but this works as an example of a hierarchical model.
code_grp_prior_mean <- root("demos_rstan", "grp_prior_mean.stan")
writeLines(readLines(code_grp_prior_mean))


#' Fit the model
#+ results='hide'
mod_grp <- cmdstan_model(stan_file = code_grp_prior_mean)
fit_grp <- mod_grp$sample(data = data_grp, seed = SEED, refresh=1000)

fit_grp$summary()
fit_grp$diagnostic_summary(diagnostics = c("divergences", "treedepth"))

#' We got a small number of divergences, so we increase
#' `adapt_delta=0.95`. Note that thisis useful only in case of small
#' number of divergences, and increasing `adapt_delta>0.99` doesn't
#' usually make sense, and in such cases there is need to investigate the
#' identifiability or parameter transformations.

#+ results='hide'
mod_grp <- cmdstan_model(stan_file = code_grp_prior_mean)
fit_grp <- mod_grp$sample(data = data_grp, seed = SEED, refresh=1000, adapt_delta=0.95)


fit_grp$summary()
fit_grp$diagnostic_summary(diagnostics = c("divergences", "treedepth"))


#' ## Unequal variance and hierarchical prior for mean and variance

code_grp_prior_mean_var <- root("demos_rstan", "grp_prior_mean_var.stan")
writeLines(readLines(code_grp_prior_mean_var))

#' Fit the model
#+ results='hide'
mod_grp <- cmdstan_model(stan_file = code_grp_prior_mean_var)
fit_grp <- mod_grp$sample(data = data_grp, seed = SEED, refresh=1000)

fit_grp$summary()
fit_grp$diagnostic_summary(diagnostics = c("divergences", "treedepth"))


#' We got a small number of divergences, so we increase
#' `adapt_delta=0.95`. Note that thisis useful only in case of small
#' number of divergences, and increasing `adapt_delta>0.99` doesn't
#' usually make sense, and in such cases there is need to investigate the
#' identifiability or parameter transformations.

#+ results='hide'
mod_grp <- cmdstan_model(stan_file = code_grp_prior_mean_var)
fit_grp <- mod_grp$sample(data = data_grp, seed = SEED, refresh=1000, adapt_delta=0.95);


fit_grp$summary()
fit_grp$diagnostic_summary(diagnostics = c("divergences", "treedepth"))


#' Plot the results
#+ warning=FALSE
temps <- fit_grp$draws(format = "df") |>
  as_tibble() |>
  select(starts_with("mu[")) |>
  setNames(c('June','July','August'))
mcmc_areas(temps) + xlab('Temperature')


#' Probabilities that June is hotter than July, June is hotter than August
#' and July is hotter than August:
#+ warning=FALSE
paste('p(TempJun > TempJul) = ', mean(temps$June > temps$July))
paste('p(TempJun > TempAug) = ', mean(temps$June > temps$August))
paste('p(TempJul > TempAug) = ', mean(temps$July > temps$August))


#' <br />
#' 
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2017-2020, Aki Vehtari, 2017 Markus Paasiniemi, licensed under BSD-3.
#' * Text &copy; 2017-2020, Aki Vehtari, licensed under CC-BY-NC 4.0.
