#' ---
#' title: "Bayesian data analysis - BRMS demos"
#' author: "Aki Vehtari"
#' date: "First version 2023-12-05. Last modified `r format(Sys.Date())`."
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
library(tibble)
library(pillar)
library(stringr)
library(brms)
options(brms.backend = "cmdstanr", mc.cores = 2)
library(posterior)
options(pillar.negative = FALSE)
library(loo)
library(priorsense)
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

#' As usual in case of generalizd linear models, (GLMs) brms defines
#' the priors on the latent model parameters. With Bernoulli the
#' default link function is logit, and thus the prior is set on
#' logit(theta). As there are no covariates logit(theta)=Intercept.
#' The brms default prior for Intercept is student_t(3, 0, 2.5), but
#' we use student_t(7, 0, 1.5) which is close to logistic
#' distribution, and thus makes the prior near-uniform for theta.
#' We can simulate from these priors to check the implied prior on theta.
#' We next compare the result to using normal(0, 1) prior on logit
#' probability. We visualize the implied priors by sampling from the priors.
data.frame(theta = plogis(ggdist::rstudent_t(n=20000, df=3, mu=0, sigma=2.5))) |>
  mcmc_hist() +
  xlim(c(0,1)) +
  labs(title='Default brms student_t(3, 0, 2.5) prior on Intercept')
data.frame(theta = plogis(ggdist::rstudent_t(n=20000, df=7, mu=0, sigma=1.5))) |>
  mcmc_hist() +
  xlim(c(0,1)) +
  labs(title='student_t(7, 0, 1.5) prior on Intercept')
#' Almost uniform prior on theta could be obtained also with normal(0,1.5)
data.frame(theta = plogis(rnorm(n=20000, mean=0, sd=1.5))) |>
  mcmc_hist() +
  xlim(c(0,1)) +
  labs(title='normal(0, 1.5) prior on Intercept')

#' Formula `y ~ 1` corresponds to a model $\mathrm{logit}(\theta) =
#\alpha\times 1 = \alpha$. `brms? denotes the $\alpha$ as `Intercept`.
#+ results='hide'
fit_bern <- brm(y ~ 1, family = bernoulli(), data = data_bern,
                prior = prior(student_t(7, 0, 1.5), class='Intercept'),
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

#' Make prior sensitivity analysis by powerscaling both prior and
#' likelihood. Focus on theta which is the quantity of interest.
theta <- draws |>
  subset_draws(variable='theta')
powerscale_sensitivity(fit_bern, prediction = \(x, ...) theta, num_args=list(digits=2)
                       )$sensitivity |>
                         filter(variable=='theta') |>
                         mutate(across(where(is.double),  ~num(.x, digits=2)))


#'
#' # Binomial model
#' 
#' Instead of sequence of 0's and 1's, we can summarize the data with
#' the number of trials and the number successes and use Binomial
#' model. The prior is specified in the 'latent space'. The actual
#' probability of success, theta = plogis(alpha), where plogis is the
#' inverse of the logistic function.
#' 
#' Binomial model with the same data and prior
data_bin <- data.frame(N = c(10), y = c(7))

#' Formula `y | trials(N) ~ 1` corresponds to a model
#' $\mathrm{logit}(\theta) = \alpha$, and the number of trials for
#' each observation is provided by `| trials(N)`
#+ results='hide'
fit_bin <- brm(y | trials(N) ~ 1, family = binomial(), data = data_bin,
               prior = prior(student_t(7, 0,1.5), class='Intercept'),
               seed = SEED, refresh = 0)

#' Check the summary of the posterior and convergence
fit_bin

#' The diagnostic indicates prior-data conflict, that is, both prior
#' and likelihood are informative. If there is true strong prior
#' information that would justify the normal(0,1) prior, then this is
#' fine, but otherwise more thinking is required (goal is not adjust
#' prior to remove diagnostic warnings withoyt thinking). In this toy
#' example, we proceed with this prior.
#' 

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
#' model for both groups and compute odds-ratio. To recreate the model
#' as two independent (separate) binomial models, we use formula `y |
#' trials(N) ~ 0 + grp2`, which corresponds to a model
#' $\mathrm{logit}(\theta) = \alpha \times 0 +
#' \beta_\mathrm{control}\times x_\mathrm{control} +
#' \beta_\mathrm{treatment}\times x_\mathrm{treatment} =
#' \beta_\mathrm{control}\times x_\mathrm{control} +
#' \beta_\mathrm{treatment}\times x_\mathrm{treatment}$, where
#' $x_\mathrm{control}$ is a vector with 1 for control and 0 for
#' treatment, and $x_\mathrm{treatemnt}$ is a vector with 1 for
#' treatemnt and 0 for control. As only of the vectors have 1, this
#' corresponds to separate models
#' $\mathrm{logit}(\theta_\mathrm{control}) = \beta_\mathrm{control}$
#' and $\mathrm{logit}(\theta_\mathrm{treatment}) =
#' \beta_\mathrm{treatment}$.  We can provide the same prior for all
#' $\beta$'s by setting the prior with `class='b'`. With prior
#' `student_t(7, 0,1.5)`, both $\beta$'s are shrunk towards 0, but
#' independently.
fit_bin2 <- brm(y | trials(N) ~ 0 + grp2, family = binomial(), data = data_bin2,
                prior = prior(student_t(7, 0,1.5), class='b'),
                seed = SEED, refresh = 0)

#' Check the summary of the posterior and convergence. brms is using
#' the first factor level `control` as the baseline and thus reports
#' the coefficient (population-level effect) for `treatment` (shown s
#' `grp2treatment`)

#' Check the summary of the posterior and convergence. With `~ 0 +
#' grp2` there is no `Intercept` and \beta_\mathrm{control} and
#' \beta_\mathrm{treatment} are presented as `grp2control` and
#' `grp2treatment`.
fit_bin2

#' Compute theta for each group and the odds-ratio. `brms` uses
#' bariable names `b_grp2control` and `b_grp2treatment` for
#' $\beta_\mathrm{control}$ and $\beta_\mathrm{treatment}$
#' respectively.
draws_bin2 <- as_draws_df(fit_bin2) |>
  mutate(theta_control = plogis(b_grp2control),
         theta_treatment = plogis(b_grp2treatment),
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

#' oddsratio 95% posterior interval
draws_bin2 |>
  subset(variable='oddsratio') |>
  summarise_draws(~quantile(.x, probs = c(0.025, 0.975)), ~mcse_quantile(.x, probs = c(0.025, 0.975)))


#' Make prior sensitivity analysis by powerscaling both prior and
#' likelihood.  Focus on oddsratio which is the quantity of
#' interest. We see that the likelihood is much more informative than
#' the prior, and we would expect to see a different posterior only
#' with a highly informative prior (possibly based on previous similar
#' experiments).
oddsratio <- draws_bin2 |>
  subset_draws(variable='oddsratio')
powerscale_sensitivity(fit_bin2, prediction = \(x, ...) oddsratio, num_args=list(digits=2)
                       )$sensitivity |>
                         filter(variable=='oddsratio') |>
                         mutate(across(where(is.double),  ~num(.x, digits=2)))

#' Above we used formula `y | trials(N) ~ 0 + grp2` to have separate
#' model for control and treatment group. An alternative model `y |
#' trials(N) ~ grp2` which is equal to `y | trials(N) ~ 1 + grp2`,
#' would correspond to a model $\mathrm{logit}(\theta) = \alpha \times
#' 1 + \beta_\mathrm{treatment}\times x_\mathrm{treatment} = \alpha +
#' \beta_\mathrm{treatment}\times x_\mathrm{treatment}. Now $\alpha$
#' models the probability of death (via logistic link) in the control
#' group and $\alpha + \beta_\mathrm{treatment}$ models the
#' probability of death (via logistic link) in the treatment
#' group. Now the models for the groups are connected. Furthermore, if
#' we set independent `student_t(7, 0, 1.5)` priors on $\alpha$ and
#' $\beta_\mathrm{treatment}$, the implied priors on
#' $\theta_\mathrm{control}$ and $\theta_\mathrm{treatment}$ are
#' different. We can verify this with a prior simulation.
#' 
data.frame(theta_control = plogis(ggdist::rstudent_t(n=20000, df=7, mu=0, sigma=1.5))) |>
  mcmc_hist() +
  xlim(c(0,1)) +
  labs(title='student_t(7, 0, 1.5) prior on Intercept') +
data.frame(theta_treatment = plogis(ggdist::rstudent_t(n=20000, df=7, mu=0, sigma=1.5))+
             plogis(ggdist::rstudent_t(n=20000, df=7, mu=0, sigma=1.5))) |>
  mcmc_hist() +
  xlim(c(0,1)) +
  labs(title='student_t(7, 0, 1.5) prior on Intercept and b_grp2treatment')

#' In this case, with relatively big treatment and control group, the
#' likelihood is informative, and the difference between using `y |
#' trials(N) ~ 0 + grp2` or `y | trials(N) ~ grp2` is negligible.
#'
#' Third option would be a hierarchical model with formula `y |
#' trials(N) ~ 1 + (1 | grp2)`, which is equivalent to `y | trials(N)
#' ~ 1 + (1 | grp2)`, and corresponds to a model
#' $\mathrm{logit}(\theta) = \alpha \times 1 +
#' \beta_\mathrm{control}\times x_\mathrm{control} +
#' \beta_\mathrm{treatment}\times x_\mathrm{treatment}$, but now the
#' prior on $\beta_\mathrm{control}$ and $\beta_\mathrm{treatment}$ is
#' $\mathrm{normal}(0, \sigma_\mathrm{grp})$. The default `brms` prior
#' for $\sigma_\mathrm{grp}$ is `student_t(3, 0, 2.5)`. Now $\alpha$
#' models the overall probablity of death (via logistic link), and
#' $\beta_\mathrm{control}$ and $\beta_\mathrm{treatment}$ model the
#' difference from that having the same prior. Prior for
#' $\beta_\mathrm{control}$ and $\beta_\mathrm{treatment}$ includes
#' unknown scale $\sigma_\mathrm{grp}$. If the there is not difference
#' between control and treatment groups, the posterior of
#' $\sigma_\mathrm{grp}$ has more mass near 0, and bigger the
#' difference between control and treatment groups are, more mass
#' there is away from 0. With just two groups, there is not much
#' information about $\sigma_\mathrm{grp}$, and unless there is a
#' informative prior on $\sigma_\mathrm{grp}$, two group hierarchical
#' model is not that useful. Hierarchical models are more useful with
#' more than two groups. In the following, we use the previously used
#' `student_t(7, 0,1.5)` prior on intercept and the default `brms`
#' prior `student_t(3, 0, 2.5)` on $\sigma_\mathrm{grp}$.
#+ results='hide'
fit_bin2 <- brm(y | trials(N) ~ 1 + (1 | grp2), family = binomial(), data = data_bin2,
                prior = prior(student_t(7, 0,1.5), class='Intercept'),
                seed = SEED, refresh = 0, control=list(adapt_delta=0.99))


#' Check the summary of the posterior and convergence. The summary
#' reports that there are Group-Level Effects: `~grp2` with 2 levels
#' (control and treatment), with `sd(Intercept)` denoting
#' $\sigma_\mathrm{grp}$. In addition, the summary lists
#' Population-Level Effects: `Intercept` ($\alpha$) as in the prevous
#' non-hierarchical models.
fit_bin2

#' We can also look at the variable names `brms` uses internally
as_draws_rvars(fit_bin2)

#' Although there is no difference, illustrate how to compute the
#' oddsratio from hierarchical model
draws_bin2 <- as_draws_df(fit_bin2)
oddsratio <- draws_bin2 |>
  mutate_variables(theta_control = plogis(b_Intercept + `r_grp2[control,Intercept]`),
                   theta_treatment = plogis(b_Intercept + `r_grp2[treatment,Intercept]`),
                   oddsratio = (theta_treatment/(1-theta_treatment))/(theta_control/(1-theta_control))) |>
  subset_draws(variable='oddsratio')
oddsratio |> mcmc_hist() +
  scale_x_continuous(breaks=seq(0.2,1.6,by=0.2))+
  geom_vline(xintercept=1, linetype='dashed')

#' Make also prior sensitivity analysis with focus on oddsratio.
powerscale_sensitivity(fit_bin2, prediction = \(x, ...) oddsratio, num_args=list(digits=2)
                       )$sensitivity |>
                         filter(variable=='oddsratio') |>
                         mutate(across(where(is.double),  ~num(.x, digits=2)))

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
#' Formula `temp ~ year` corresponds to model $\mathrm{temp} ~
#' \mathrm{normal}(\alpha + \beta \times \mathrm{temp}, \sigma).  The
#' model could also be defined as `temp ~ 1 + year` which explicitly
#' shows the intercept ($\alpha$) part. Using the variable names
#' `brms` uses the model can be written also as `temp ~
#' normal(b_Intercept*1 + b_year*year, sigma)`. We start with the
#' default priors to see some tricks that `brms` does behind the
#' curtain.
fit_lin <- brm(temp ~ year, data = data_lin, family = gaussian(),
               seed = SEED, refresh = 0)

#' Check the summary of the posterior and convergence.
fit_lin

#' Convergence diagnostics look good. We see that posterior mean of
#' `Intercept` is -34.7, which may sound strange, but that is the
#' intercept at year 0, that is, very far from the data range, and
#' thus doesn't have meaningful interpretation directly. The posterior
#' mean of `year` coefficient is 0.02, that is, we estimate that the
#' summer temperature is increasing 0.02°C per year (which would make
#' 1°C in 50 years).
#' 
#' We can check $R^2$ which corresponds to the proporion of variance
#' explained by the model. The linear model explains 0.16=16% of the
#' total data variance.
bayes_R2(fit_lin) |> round(2)

#' We can check the all the priors used. 
prior_summary(fit_lin)

#' We see that `class=b` and `coef=year` have `flat`, that is,
#' improper uniform prior, `Intercept` has `student_t(3, 9.5, 2.5)`,
#' and `sigma` has `student_t(3, 0, 2.5)` prior.  In general it is
#' good to use proper priors, but sometimes flat priors are fine and
#' produce proper posterior (like in this case). Important part here
#' is that by default, `brms` sets the prior on Intercept after
#' centering the covariate values (design matrix). In this case,
#' `brms` uses `temp - mean(temp) = temp - 1987` instead of original
#' years. This in general improves the sampling efficiency. As the
#' `Intercept` is now defined at the middle of the data, the default
#' `Intercept` prior is centered on median of the target (here target
#' is `year`). If we would like to set informative priors, we need to
#' set the informative prior on `Intercept` given the centered
#' covariate values. We can turn of the centering by setting argument
#' `center=FALSE`, and we can set the prior on original intercept by
#' using a formula `temp ~ 0 + Intercept + year`. In this case, we are
#' happy with the default prior for the intercept. In this specific
#' casse, the flat prior on coefficient is also fine, but we add an
#' weakly informative prior just for the illustration. Let's assume we
#' expect the temperature to change less than 1°C in 10 years. With
#' `student_t(3, 0, 0.03)` about 95% prior mass has less than 0.1°C
#' change in year, and with low degrees of freedom (3) we have thick
#' tails making the likelihood dominate in case of prior-data
#' conflict. In real life, we do have much more information about the
#' temperature change, and naturally a hierarchical spatio-temporal
#' model with all temperature measurement locations would be even
#' better.
fit_lin <- brm(temp ~ year, data = data_lin, family = gaussian(),
               prior = prior(student_t(3, 0, 0.03), class='b'),
               seed = SEED, refresh = 0)

#' Check the summary of the posterior and convergence
fit_lin

#' Make prior sensitivity analysis by powerscaling both prior and likelihood.
powerscale_sensitivity(fit_lin)$sensitivity |>
                                mutate(across(where(is.double),  ~num(.x, digits=2)))

#' Our weakly informative proper prior has negligible sensitivity, and
#' the likelihood is informative.

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
#' All posterior draws have `b_year>0`, the probability gets rounded
#' to 1, and MCSE is not available as the obserevd posterior variance
#' is 0.
#' 

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

#' # Linear Student's $t$ model
#' 
#' The temperatures used in the above analyses are averages over three
#' months, which makes it more likely that they are normally
#' distributed, but there can be extreme events in the feather and we
#' can check whether more robust Student's $t$ observation model would
#' give different results.
#' 
#+  results='hide'
fit_lin_t <- brm(temp ~ year, data = data_lin, family = student(),
                 prior = prior(student_t(3, 0, 0.03), class='b'),
                 seed = SEED, refresh = 0)

#' Check the summary of the posterior and convergence. The b_year
#' posterior looks similar as before and the posterior for degrees of
#' freedom `nu` has most of the posterior mass for quite large values
#' indicating there is no strong support for thick tailed variation in
#' average summer temperatures.
fit_lin_t

#' # Pareto-smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO)
#' 
#' We can use leave-one-out cross-validation to compare the expected predictive performance.
#' 
#' LOO comparison shows normal and Student's $t$ model have similar performance.
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
                 prior = prior(student_t(3, 0, 0.03), class='b'),
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

#' Make prior sensitivity analysis by powerscaling both prior and likelihood.
powerscale_sensitivity(fit_lin_h)$sensitivity |>
                                mutate(across(where(is.double),  ~num(.x, digits=2)))


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
fit_spline_h <- brm(bf(temp ~ s(year),
                     sigma ~ s(year)),
                  data = data_lin, family = gaussian(),
                  seed = SEED, refresh = 0)

#' We get warnings about divergences, and try rerunning with higher
#' adapt_delta, which leads to using smaller step sizes. Often
#' `adapt_delta=0.999` leads to very slow sampling, but with this
#' small data, this is not an issue.
fit_spline_h <- update(fit_spline_h, control = list(adapt_delta=0.999))

#' Check the summary of the posterior and convergence. We're not
#' anymore able to make interpretation of the temperature increase
#' based on this summary. For splines, we see prior scales `sds` for
#' the spline coefficients.
fit_spline_h

#' We can still plot posterior predictive distribution at each year
#' until 2030 `add_predicted_draws()` takes the years from the data
#' and uses `fit_lin_h` to make the predictions.
data_lin |>
  add_row(year=2023:2030) |>
  add_predicted_draws(fit_spline_h) |>
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

#' And we can use leave-one-out cross-validation to compare the
#' expected predictive performance.
#' 
#' LOO comparison shows homoskedastic normal linear and
#' heteroskedastic normal spline models have similar
#' performances. There are not enough observations to make clear
#' difference between the models.
loo_compare(loo(fit_lin), loo(fit_spline_h))

#' For spline and other non-parametric models, we can use predictive
#' estimates and predictions to get interpretable quantities. Let's
#' examine the difference of estimated average temperature in years
#' 1952 and 2022.
temp_diff <- posterior_epred(fit_spline_h, newdata=filter(data_lin,year==1952|year==2022)) |>
  rvar() |>
  diff() |>
  as_draws_df() |>
  set_variables('temp_diff')

temp_diff <- data_lin |>
  filter(year==1952|year==2022) |>
  add_epred_draws(fit_spline_h) |>
  pivot_wider(id_cols=.draw, names_from = year, values_from = .epred) |>
  mutate(temp_diff = `2022`-`1952`,
         .chain = (.draw - 1) %/% 1000 + 1,
         .iteration = (.draw - 1) %% 1000 + 1) |>
  as_draws_df() |>
  subset_draws(variable='temp_diff')

#' Posterior distribution for average summer temperature increase from 1952 to 2022
temp_diff |>
  mcmc_hist()

#' 95% posterior interval for average summer temperature increase from 1952 to 2022
temp_diff |>
  summarise_draws(~quantile(.x, probs = c(0.025, 0.975)),
                  ~mcse_quantile(.x, probs = c(0.025, 0.975)),
                  .num_args = list(digits = 2, notation = "dec"))

#' Make prior sensitivity analysis by powerscaling both prior and
#' likelihood with focus on average summer temperature increase from
#' 1952 to 2022.
powerscale_sensitivity(fit_spline_h, prediction = \(x, ...) temp_diff, num_args=list(digits=2)
                       )$sensitivity |>
                         filter(variable=='temp_diff') |>
                         mutate(across(where(is.double),  ~num(.x, digits=2)))

#' Probability that the average summer temperature has increased from
#' 1952 to 2022 is 99.5%.
temp_diff |>
  mutate(I_temp_diff_gt_0 = temp_diff>0,
         temp_diff = NULL) |>
  subset_draws(variable='I_temp_diff_gt_0') |>
  summarise_draws(mean, mcse_mean)


#' 
#' # Comparison of k groups with hierarchical normal models
#'
#' Load factory data, which contain 5 quality measurements for each of
#' 6 machines. We're interested in analysing are the quality differences
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
#+ results='hide'
fit_pooled <- brm(quality ~ 1, data = factory, refresh=0)

#' Check the summary of the posterior and convergence.
fit_pooled

#' ## Separate model
#'
#' As comparison make also seprate model. To make it completely
#' separate we need to have different sigma for each machine, too.
#+ results='hide'
fit_separate <- brm(bf(quality ~ 0 + machine,
                       sigma ~ 0 + machine),
                    data = factory, refresh=0)

#' Check the summary of the posterior and convergence.
fit_separate

#' # Common variance hierarchical model (ANOVA)
#+ results='hide'
fit_hier <- brm(quality ~ 1 + (1 | machine),
                data = factory, refresh = 0)

#' Check the summary of the posterior and convergence.
fit_hier

#' LOO comparison shows the hierarchical model is the best. The
#' differences are small as the number of observations is small and
#' there is a considerable prediction (aleatoric) uncertainty.
loo_compare(loo(fit_pooled), loo(fit_separate), loo(fit_hier))

#' Different model posterior distributions for the mean
#' quality. Pooled model ignores the varition between
#' machines. Separate model doesn't take benefit from the similariy of
#' the machines and has higher uncertainty.
ph <- fit_hier |>
  spread_rvars(b_Intercept, r_machine[machine,]) |>
  mutate(machine_mean = b_Intercept + r_machine) |>
  ggplot(aes(xdist=machine_mean, y=machine)) +
  stat_halfeye() +
  scale_y_continuous(breaks=1:6) +
  labs(x='Quality', y='Machine', title='Hierarchical')

ps <- fit_separate |>
  as_draws_df() |>
  subset_draws(variable='b_machine', regex=TRUE) |>
  set_variables(paste0('b_machine[', 1:6, ']')) |>
  as_draws_rvars() |>
  spread_rvars(b_machine[machine]) |>
  mutate(machine_mean = b_machine) |>
  ggplot(aes(xdist=machine_mean, y=machine)) +
  stat_halfeye() +
  scale_y_continuous(breaks=1:6) +
  labs(x='Quality', y='Machine', title='Separate')

pp <- fit_pooled |>
  spread_rvars(b_Intercept) |>
  mutate(machine_mean = b_Intercept) |>
  ggplot(aes(xdist=machine_mean, y=0)) +
  stat_halfeye() +
  scale_y_continuous(breaks=NULL) +
  labs(x='Quality', y='All machines', title='Pooled')

(pp / ps / ph) * xlim(c(50,140))

#' Make prior sensitivity analysis by powerscaling both prior and
#' likelihood with focus on mean quality of each machine. We see no
#' prior sensitivity.
machine_mean <- fit_hier |>
  as_draws_df() |>
  mutate(across(matches('r_machine'), ~ .x - b_Intercept)) |>
  subset_draws(variable='r_machine', regex=TRUE) |>
  set_variables(paste0('machine_mean[', 1:6, ']'))
powerscale_sensitivity(fit_hier, prediction = \(x, ...) machine_mean, num_args=list(digits=2)
                       )$sensitivity |>
                         filter(str_detect(variable,'machine_mean')) |>
                         mutate(across(where(is.double),  ~num(.x, digits=2)))

#' 
#' # Hierarchical binomial model
#'
#' [Sorafenib Toxicity Dataset in `metadat` R package](https://wviechtb.github.io/metadat/reference/dat.ursino2021.html)
#' includes results frm 13 studies investigating the occurrence of
#' dose limiting toxicities (DLTs) at different doses of Sorafenib.
#'
#' Load data
load(url('https://github.com/wviechtb/metadat/raw/master/data/dat.ursino2021.rda'))
head(dat.ursino2021)

#' Pooled model assumes all studies have the same dose effect
#' (reminder: `~ dose` is equivalent to `~ 1 + dose`)
#+ results='hide'
fit_pooled <- brm(events | trials(total) ~ dose,
                  prior = c(prior(student_t(7, 0, 1.5), class='Intercept'),
                            prior(normal(0, 1), class='b')),
                  family=binomial(), data=dat.ursino2021)

#' Check the summary of the posterior and convergence
fit_pooled

#' Dose coefficient seems to be very small. Looking at the posterior,
#' we see that it is positive with high probability. 
fit_pooled |>
  as_draws() |>
  subset_draws(variable='b_dose') |>
  summarise_draws(~quantile(.x, probs = c(0.025, 0.975)), ~mcse_quantile(.x, probs = c(0.025, 0.975)))

#' The dose was reported in mg, and most values are in hundreds. It is
#' often sensible to switch to a scale in which the range of values is
#' closer to unit range. In this case it is natural to use g instead
#' of mg.
dat.ursino2021 <- dat.ursino2021 |>
  mutate(doseg = dose/100)

#' Fit the pooled model again uing `doseg`
#+ results='hide'
fit_pooled <- brm(events | trials(total) ~ doseg,
                  prior = c(prior(student_t(7, 0, 1.5), class='Intercept'),
                            prior(normal(0, 1), class='b')),
                  family=binomial(), data=dat.ursino2021)

#' Check the summary of the posterior and convergence.
fit_pooled

#' Now it is easier to interpret the presented values. 

#' Separate model assumes all studies have different dose effect.
#' It would be a bit complicated to set a different prior on study specific
#' intercepts and other coefficients, so we use the ame prior for all.
#+ results='hide'
fit_separate <- brm(events | trials(total) ~ 0 + study + doseg:study,
                    prior=prior(student_t(7, 0, 1.5), class='b'),
                    family=binomial(), data=dat.ursino2021)

#' Check the summary of the posterior and convergence.
fit_separate

#' Hierarchical model assumes common mean effect and variation around with normal population prior
#' (reminder: `~ doseg + (doseg | study)` is equivalent to `~ 1 + doseg + (1 + doseg | study)`)
#+ results='hide'
fit_hier <- brm(events | trials(total) ~ doseg + (doseg | study),
                    prior=c(prior(student_t(7, 0, 1.5), class='Intercept'),
                            prior(normal(0, 1), class='b')),
                family=binomial(), data=dat.ursino2021)

#' We seem some divergences and repeat with higher adapt_delta
#+ results='hide'
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
mcmc_areas(as_draws_df(fit_hier), regex_pars='r_study\\[.*doseg')
#'
#' The population level coefficient for the dose is clearly larger than 0
mcmc_areas(as_draws_df(fit_hier), regex_pars='b_doseg') +
  geom_vline(xintercept=0, linetype='dashed') +
  xlim(c(0,1))

#' Make prior sensitivity analysis by powerscaling both prior and likelihood.
powerscale_sensitivity(fit_hier, variable='b_doseg'
                       )$sensitivity |>
                         mutate(across(where(is.double),  ~num(.x, digits=2)))

#' The posterior for the probability of event given certain dose and a new study.
data.frame(study='new',
           doseg=seq(0.1,1,by=0.1),
           total=1) |>
  add_linpred_draws(fit_hier, transform=TRUE, allow_new_levels=TRUE) |>
  ggplot(aes(x=doseg, y=.linpred)) +
  stat_lineribbon(.width = c(.95), alpha = 1/2, color=brewer.pal(5, "Blues")[[5]]) +
  scale_fill_brewer()+
  labs(x= "Dose (g)", y = 'Probability of event') +
  theme(legend.position="none") +
  geom_hline(yintercept=0) +
  scale_x_continuous(breaks=seq(0.1,1,by=0.1))

#' If plot individual posterior draws, we see that there is a lot of
#' uncertainty about the overall probability (explained by the
#' variation in Intercept in different studies), but less uncertainty
#' about the slope.
data.frame(study='new',
           doseg=seq(0.1,1,by=0.1),
           total=1) |>
  add_linpred_draws(fit_hier, transform=TRUE, allow_new_levels=TRUE, ndraws=100) |>
  ggplot(aes(x=doseg, y=.linpred)) +
  geom_line(aes(group=.draw), alpha = 1/2, color = brewer.pal(5, "Blues")[[3]])+
  scale_fill_brewer()+
  labs(x= "Dose (g)", y = 'Probability of event') +
  theme(legend.position="none") +
  geom_hline(yintercept=0) +
  scale_x_continuous(breaks=seq(0.1,1,by=0.1))

#' Posterior predictive checking showing the observed and predicted number of events.
pp_check(fit_hier, type = "ribbon_grouped", group="study")

#' <br />
#' 
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2017-2024, Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2017-2024, Aki Vehtari, licensed under CC-BY-NC 4.0.
