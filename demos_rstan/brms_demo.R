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
library(brms)
options(brms.backend = "cmdstanr", mc.cores = 2)
library(posterior)
options(pillar.negative = FALSE)
library(loo)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(tidybayes)
library(ggdist)
library(patchwork)
library(RColorBrewer)
SEED <- 48927 # set random seed for reproducability

#' # Introduction
#' 
#' This notebook contains several examples of how to use [Stan](https://mc-stan.org) in R with [__brms__](https://paul-buerkner.github.io/brms/). This notebook assumes basic knowledge of Bayesian inference and MCMC. The examples are related to [Bayesian data analysis course](https://avehtari.github.io/BDA_course_Aalto/).
#' 
#' # Bernoulli model
#' 
#' Toy data with sequence of failures (0) and successes (1). We would
#' like to learn about the unknown probability of success.
data_bern <- data.frame(y = c(1, 1, 1, 0, 1, 1, 1, 0, 1, 0))


#' brms uses by default student_t(3, 0, 2.5), bu we can assign uniform
#' prior (beta(1,1)). 
#+ results='hide'
fit_bern <- brm(y ~ 1, family = bernoulli(), data = data_bern,
                prior = prior("", class='Intercept'),
                seed = SEED, refresh = 0)

#' Check the summary of the posterior and convergence
fit_bern

#' Extract the posterior draws
draws <- as_draws_df(fit_bern)
#' We can get summary information using summarise_draws()
draws |>
  subset_draws(variable='b_Intercept') |>
  summarise_draws()

#' We can compute the probability of success by using plogis which is
#' equal to inverse-logit function
draws <- draws |>
  mutate_variables(theta=plogis(b_Intercept))
#' Summary of theta by using summarise_draws()
draws |>
  subset_draws(variable='theta') |>
  summarise_draws()

#' Histogram of theta
mcmc_hist(draws, pars='theta') +
  xlab('theta') +
  xlim(c(0,1))

#' We next compare the result to using normal(0, 1) prior on logit
#' probability. Visualize the prior by drawing samples from it
prior_mean <- 0
prior_sd <- 1
prior_intercept <- normal(location = prior_mean, scale = prior_sd)
prior_samples <- data.frame(
                 theta = plogis(rnorm(20000, prior_mean, prior_sd)))
mcmc_hist(prior_samples) +
  xlim(c(0,1))

#+ results='hide'
fit_bern <- brm(y ~ 1, family = bernoulli(), data = data_bern,
                prior = prior(normal(0, 1), class='Intercept'),
                seed = SEED, refresh = 0)

#' Check the summary of the posterior and convergence
fit_bern

#' We can examine the latent parameter
draws <- as_draws_df(fit_bern)
draws |>
  subset_draws(variable='b_Intercept') |>
  summarise_draws()

#' We can compute the probability of success by using plogis which is
#' equal to inverse-logit function
draws <- draws |>
  mutate_variables(theta=plogis(b_Intercept))
#' Summary of theta by using summarise_draws()
draws |>
  subset_draws(variable='theta') |>
  summarise_draws()

#' Histogram of theta
mcmc_hist(draws, pars='theta') +
  xlab('theta') +
  xlim(c(0,1))

#' As the number of observations is small, there is small change in
#' the posterior mean when the prior is changed. You can experiment
#' with different priors and varying the number of observations.
#' 
#' 
#' # Binomial model
#' 
#' Instead of sequence of 0's and 1's, we can summarize the data with
#' the number of trials and the number successes and use Binomial
#' model. The prior is specified in the 'latent space'. The actual
#' probability of success, theta = plogis(alpha), where plogis is the
#' inverse of the logistic function.
#' 
#' Visualize the prior by drawing samples from it
prior_mean <- 0
prior_sd <- 1
prior_samples <- data.frame(
                 theta = plogis(rnorm(20000, prior_mean, prior_sd)))
mcmc_hist(prior_samples)

#' Binomial model with the same data
data_bin <- data.frame(N = c(10), y = c(7))
#+ results='hide'
fit_bin <- brm(y | trials(N) ~ 1, family = binomial(), data = data_bin,
               prior = prior(normal(0,1), class='Intercept'),
               seed = SEED, refresh = 0)

#' Check the summary of the posterior and convergence
fit_bin

#' Extract the posterior draws
draws <- as_draws_df(fit_bin)
#' We can get summary information using summarise_draws()
draws |>
  subset_draws(variable='b_Intercept') |>
  summarise_draws()

#' We can compute the probability of success by using plogis which is
#' equal to inverse-logit function
draws <- draws |>
  mutate_variables(theta=plogis(b_Intercept))
#' Summary of theta by using summarise_draws()
draws |>
  subset_draws(variable='theta') |>
  summarise_draws()

#' Histogram of theta
mcmc_hist(draws, pars='theta') +
  xlab('theta') +
  xlim(c(0,1))

#' Re-run the model with a new data dataset without recompiling
data_bin <- data.frame(N = c(5), y = c(4))
fit_bin <- update(fit_bin, newdata = data_bin)

#' Check the summary of the posterior and convergence
fit_bin

#' Extract the posterior draws
draws <- as_draws_df(fit_bin)
#' We can get summary information using summarise_draws()
draws |>
  subset_draws(variable='b_Intercept') |>
  summarise_draws()

#' We can compute the probability of success by using plogis which is
#' equal to inverse-logit function
draws <- draws |>
  mutate_variables(theta=plogis(b_Intercept))
#' Summary of theta by using summarise_draws()
draws |>
  subset_draws(variable='theta') |>
  summarise_draws()

#' Histogram of theta
mcmc_hist(draws, pars='theta') +
  xlab('theta') +
  xlim(c(0,1))

#' # Comparison of two groups with Binomial 
#' 
#' An experiment was performed to estimate the effect of beta-blockers
#' on mortality of cardiac patients. A group of patients were randomly
#' assigned to treatment and control groups:
#' 
#' - out of 674 patients receiving the control, 39 died
#' - out of 680 receiving the treatment, 22 died
#' 
#' Data, where `grp2` is an indicator variable defined as a factor
#' type, which is useful for categorical variables.
data_bin2 <- data.frame(N = c(674, 680), y = c(39,22), grp2 = factor(c('control','treatment')))

#' To analyse whether the treatment is useful, we can use Binomial
#' model for both groups and compute odds-ratio.
fit_bin2 <- brm(y | trials(N) ~ grp2, family = binomial(), data = data_bin2,
                prior = prior(normal(0,1), class='Intercept'),
                seed = SEED, refresh = 0)

#' Check the summary of the posterior and convergence. brms is using
#' the first factor level `control` as the baseline and thus reports
#' the coefficient (population-level effect) for `treatment` (shown s
#' `grp2treatment`)
fit_bin2

#' Compute theta for each group and the odds-ratio
draws_bin2 <- as_draws_df(fit_bin2) |>
  mutate(theta_control = plogis(b_Intercept),
         theta_treatment = plogis(b_Intercept + b_grp2treatment),
         oddsratio = (theta_treatment/(1-theta_treatment))/(theta_control/(1-theta_control)))
#' Plot oddsratio
mcmc_hist(draws_bin2, pars='oddsratio') +
  scale_x_continuous(breaks=seq(0.2,1.6,by=0.2))+
  geom_vline(xintercept=1, linetype='dashed')

#' Probability that the oddsratio<1
draws_bin2 |>
  mutate(poddsratio = oddsratio<1) |>
  subset(variable='poddsratio') |>
  summarise_draws(mean, mcse_mean)

#' oddratio 95% posterior interval
draws_bin2 |>
  subset(variable='oddsratio') |>
  summarise_draws(~quantile(.x, probs = c(0.025, 0.975)), ~mcse_quantile(.x, probs = c(0.025, 0.975)))

#' # Linear Gaussian model
#' 
#' Use the Kilpisjärvi summer month temperatures 1952--2022 data from `aaltobda` package
load(url('https://github.com/avehtari/BDA_course_Aalto/raw/master/rpackage/data/kilpisjarvi2022.rda'))
data_lin <- data.frame(year = kilpisjarvi2022$year,
                       temp = kilpisjarvi2022$temp.summer)

#' Plot the data
data_lin |>
  ggplot(aes(year, temp)) +
  geom_point(color=2) +
  labs(x= "Year", y = 'Summer temp. @Kilpisjärvi') +
  guides(linetype = "none")

#' To analyse has there been change in the average summer month
#' temperature we use a linear model with Gaussian model for the
#' unexplained variation. By default brms uses uniform prior for the
#' coefficients.
#' 
#' `temp ~ year` means temp depends on the intercept and `temp`.
#' The model could also be defined as `temp ~ 1 + year` which explicitly shows the
#' intercept part. The corresponding regression model is
#' temp ~ normal(b_Intercept*1 + b_year*year, sigma)
fit_lin <- brm(temp ~ year, data = data_lin, family = gaussian(),
               seed = SEED, refresh = 0)

#' We can check the all the priors used. In general it is good to use
#' proper priors, but sometimes flat priors are fine and produce
#' proper posterior.
prior_summary(fit_lin)

#' Check the summary of the posterior and convergence
fit_lin

#' Extract the posterior draws and check the summaries
draws_lin <- as_draws_df(fit_lin) 
draws_lin |> summarise_draws()
#' If one of the columns is hidden we can force printing all columns
draws_lin |> summarise_draws() |> print(width=Inf)

#' Histogram of b_year
draws_lin |>
  mcmc_hist(pars='b_year') +
  xlab('Average temperature increase per year')

#' Probability that the coefficient b_year > 0 and the corresponding MCSE
draws_lin |>
  mutate(I_b_year_gt_0 = b_year>0) |>
  subset_draws(variable='I_b_year_gt_0') |>
  summarise_draws(mean, mcse_mean)

#' 95% posterior interval for temperature increase per 100 years
draws_lin |>
  mutate(b_year_100 = b_year*100) |>
  subset_draws(variable='b_year_100') |>
  summarise_draws(~quantile(.x, probs = c(0.025, 0.975)),
                  ~mcse_quantile(.x, probs = c(0.025, 0.975)),
                  .num_args = list(digits = 2, notation = "dec"))

#' Plot posterior draws of the linear function values at each year.
#' `add_linpred_draws()` takes the years from the data and uses `fit_lin` to make
#' the predictions.
data_lin |>
  add_linpred_draws(fit_lin) |>
  # plot data
  ggplot(aes(x=year, y=temp)) +
  geom_point(color=2) +
  # plot lineribbon for the linear model
  stat_lineribbon(aes(y = .linpred), .width = c(.95), alpha = 1/2, color=brewer.pal(5, "Blues")[[5]]) +
  # decoration
  scale_fill_brewer()+
  labs(x= "Year", y = 'Summer temp. @Kilpisjärvi') +
  theme(legend.position="none")+
  scale_x_continuous(breaks=seq(1950,2020,by=10))

#' Alternativelly plot a spaghetti plot for 100 draws
data_lin |>
  add_linpred_draws(fit_lin, ndraws=100) |>
  # plot data
  ggplot(aes(x=year, y=temp)) +
  geom_point(color=2) +
  # plot a line for each posterior draw
  geom_line(aes(y=.linpred, group=.draw), alpha = 1/2, color = brewer.pal(5, "Blues")[[3]])+
  # decoration
  scale_fill_brewer()+
  labs(x= "Year", y = 'Summer temp. @Kilpisjärvi') +
  theme(legend.position="none")+
  scale_x_continuous(breaks=seq(1950,2020,by=10))

#' Plot posterior predictive distribution at each year until 2030
#' `add_predicted_draws()` takes the years from the data and uses
#' `fit_lin` to make the predictions.
data_lin |>
  add_row(year=2023:2030) |>
  add_predicted_draws(fit_lin) |>
  # plot data
  ggplot(aes(x=year, y=temp)) +
  geom_point(color=2) +
  # plot lineribbon for the linear model
  stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 1/2, color=brewer.pal(5, "Blues")[[5]]) +
  # decoration
  scale_fill_brewer()+
  labs(x= "Year", y = 'Summer temp. @Kilpisjärvi') +
  theme(legend.position="none")+
  scale_x_continuous(breaks=seq(1950,2030,by=10))

#' # Linear Student's t model
#' 
#' The temperatures used in the above analyses are averages over three
#' months, which makes it more likely that they are normally
#' distributed, but there can be extreme events in the feather and we
#' can check whether more robust Student's t observation model woul
#' give different results.
#' 
#+  results='hide'
fit_lin_t <- brm(temp ~ year, data = data_lin, family = student(),
                 seed = SEED, refresh = 0)

#' Check the summary of the posterior and convergence. The b_year
#' posterior looks similar as before and the posterior for degrees of
#' freedom `nu` has most of the posterior mas for quite large values
#' indicating there is no strong support for thick tailed variation in
#' temperature.
fit_lin_t

#' # Pareto-smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO)
#' 
#' We can use leave-one-out cross-validation to compare the expected predictive performance.
#' 
#' LOO comparison shows normal and Student's t model have similar performance.
loo_compare(loo(fit_lin), loo(fit_lin_t))

#' # Heteroskedastic linear model
#' 
#' Heteroskedasticity assumes that the variation around the linear
#' mean can also vary. We can allow sigma to depend on year, too.
#' Although the additional component is written as `sigma ~ year`, the
#' log link function is used and the model is for log(sigma). `bf()` allows
#' listing several formulas.
#' 
#+  results='hide'
fit_lin_h <- brm(bf(temp ~ year,
                    sigma ~ year),
                 data = data_lin, family = gaussian(),
                 seed = SEED, refresh = 0)

#' Check the summary of the posterior and convergence. The b_year
#' posterior looks similar as before. The posterior for sigma_year
#' looks like having mosst of the ma for negative values, indicating
#' decrease in temperature variation around the mean.
fit_lin_h

#' Histogram of b_year and b_sigma_year
as_draws_df(fit_lin_h) |>
  mcmc_areas(pars=c('b_year', 'b_sigma_year'))

#' As log(x) is almost linear when x is close to zero, we can see that the
#' sigma is decreasing about 1% per year (95% interval from 0% to 2%).
#' 

#' Plot posterior predictive distribution at each year until 2030
#' `add_predicted_draws()` takes the years from the data and uses
#' `fit_lin_h` to make the predictions.
data_lin |>
  add_row(year=2023:2030) |>
  add_predicted_draws(fit_lin_h) |>
  # plot data
  ggplot(aes(x=year, y=temp)) +
  geom_point(color=2) +
  # plot lineribbon for the linear model
  stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 1/2, color=brewer.pal(5, "Blues")[[5]]) +
  # decoration
  scale_fill_brewer()+
  labs(x= "Year", y = 'Summer temp. @Kilpisjärvi') +
  theme(legend.position="none")+
  scale_x_continuous(breaks=seq(1950,2030,by=10))

#' We can use leave-one-out cross-validation to compare the expected predictive performance.
#' 
#' LOO comparison shows homoskedastic normal and heteroskedastic
#' normal models have similar performances.
loo_compare(loo(fit_lin), loo(fit_lin_h))

#' # Heteroskedastic non-linear model
#' 
#' We can test the linearity assumption by using non-linear spline
#' functions, by uing `s(year)` terms. Sampling is slower as the
#' posterior gets more complex.
#' 
#+  results='hide'
fit_lin_hs <- brm(bf(temp ~ s(year),
                     sigma ~ s(year)),
                  data = data_lin, family = gaussian(),
                  seed = SEED, refresh = 0)

#' We get warnings about divergences, and try rerunning with higher
#' adapt_delta, which leads to using smaller step sizes. Often
#' `adapt_delta=0.999` leads to very slow sampling, but with this
#' small data, this is not an issue.
fit_lin_hs <- update(fit_lin_hs, control = list(adapt_delta=0.999))

#' Check the summary of the posterior and convergence. The b_year
#' posterior looks similar as before. The posterior for sigma_year
#' looks like having mosst of the ma for negative values, indicating
#' decrease in temperature variation around the mean.
fit_lin_hs

#' Plot posterior predictive distribution at each year until 2030
#' `add_predicted_draws()` takes the years from the data and uses
#' `fit_lin_h` to make the predictions.
data_lin |>
  add_row(year=2023:2030) |>
  add_predicted_draws(fit_lin_hs) |>
  # plot data
  ggplot(aes(x=year, y=temp)) +
  geom_point(color=2) +
  # plot lineribbon for the linear model
  stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 1/2, color=brewer.pal(5, "Blues")[[5]]) +
  # decoration
  scale_fill_brewer()+
  labs(x= "Year", y = 'Summer temp. @Kilpisjärvi') +
  theme(legend.position="none")+
  scale_x_continuous(breaks=seq(1950,2030,by=10))

#' We can use leave-one-out cross-validation to compare the expected predictive performance.
#' 
#' LOO comparison shows homoskedastic normal linear and
#' heteroskedastic normal spline models have similar
#' performances. There are not enough observations to make clear
#' difference between the models.
loo_compare(loo(fit_lin), loo(fit_lin_hs))

#' 
#' # Comparison of k groups with hierarchical normal models
#'
#' Load factory data, which contain 5 quality measurements for each of
#' 6 machines. We're interested in analying are the quality differences
#' between the machines.
factory <- read.table(url('https://raw.githubusercontent.com/avehtari/BDA_course_Aalto/master/rpackage/data-raw/factory.txt'))
colnames(factory) <- 1:6
factory
#' We pivot the data to long format
factory <- factory |>
  pivot_longer(cols = everything(),
               names_to = 'machine',
               values_to = 'quality')
factory

#' ## Pooled model
#'
#' As comparison make also pooled model
fit_pooled <- brm(quality ~ 1, data = factory, refresh=0)

#' ## Separate model
#'
#' As comparison make also seprate model. To make it completely
#' separate we need to have different sigma for each machine, too.
fit_separate <- brm(bf(quality ~ machine,
                       sigma ~ machine),
                    data = factory, refresh=0)

#' # Common variance hierarchical model (ANOVA)
#' 
fit_hier <- brm(quality ~ 1 + (1 | machine),
                data = factory, refresh = 0)

#' Check the summary of the posterior and convergence.
fit_hier

#' LOO comparison shows the hierarchical model is the best
loo_compare(loo(fit_pooled), loo(fit_separate), loo(fit_hier))

#' Distributions of quality differences from the mean quality
mcmc_areas(as_draws_df(fit_hier), regex_pars='r_machine')

#' Posterior predictive distributions for 6 old and 1 new machines
posterior_predict(fit_hier, newdata=data.frame(machine=1:7, quality=rep(NA,7)),
                  allow_new_levels=TRUE) |>
  as_draws_df() |>
  mcmc_areas()

#' 
#' # Hierarchical model
#'
#' 
#'

load(url('https://github.com/wviechtb/metadat/raw/master/data/dat.ursino2021.rda'))
head(dat.ursino2021)

#' Pooled model assumes all studies have the same dose effect
fit_pooled <- brm(events | trials(total) ~ dose,
                  family=binomial(), data=dat.ursino2021)

#' Separate model assumes all studies have different dose effect
fit_separate <- brm(events | trials(total) ~ dose:study,
                    family=binomial(), data=dat.ursino2021)
fit_separate <- update(fit_separate, control=list(init=0.1))

#' Hierarchical model assumes common mean effect and variation round with normal population prior
fit_hier <- brm(events | trials(total) ~ dose + (dose | study),
                family=binomial(), data=dat.ursino2021)
fit_hier <- update(fit_hier, control=list(adapt_delta=0.99))

#' LOO-CV comparison
loo_compare(loo(fit_pooled), loo(fit_separate), loo(fit_hier))

#' We get warnings about Pareto k's > 0.7 in PSIS-LOO for separate
#' model, but as in that case the LOO-CV estimate is usually
#' overoptimistic and the separate model is the worst, there is no
#' need to use more accurate computation.
#'
#' Hierarchical model has better elpd than the pooled, but difference
#' is negligible. However, when we look at the study specific
#' parameters, we see that the Miller study has higher intercept (more
#' events).
mcmc_areas(as_draws_df(fit_hier), regex_pars='r_study\\[.*Intercept')
#'
#' There are no differences in slopes.
mcmc_areas(as_draws_df(fit_hier), regex_pars='r_study\\[.*dose')
#'
#' The coefficient for the dose is clearly way from 0
mcmc_areas(as_draws_df(fit_hier), regex_pars='b_dose') +
  geom_vline(xintercept=0, linetype='dashed') +
  xlim(c(0,0.01))

#' The posterior for the probability of event given certain dose and a new study
data.frame(study='new',
           dose=seq(100,1000,by=100),
           total=1) |>
  add_linpred_draws(fit_hier, transform=TRUE, allow_new_levels=TRUE) |>
  ggplot(aes(x=dose, y=.linpred)) +
  stat_lineribbon(.width = c(.95), alpha = 1/2, color=brewer.pal(5, "Blues")[[5]]) +
  scale_fill_brewer()+
  labs(x= "Dose", y = 'Probability of event') +
  theme(legend.position="none") +
  geom_hline(yintercept=0) +
  scale_x_continuous(breaks=seq(100,1000,by=100))

#' Posterior predictive checking
pp_check(fit_hier, type = "ribbon_grouped", group="study")

#' <br />
#' 
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2017-2023, Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2017-2023, Aki Vehtari, licensed under CC-BY-NC 4.0.
