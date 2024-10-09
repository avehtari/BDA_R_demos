#' ---
#' title: "Bayesian data analysis - BRMS demos"
#' author: "Aki Vehtari"
#' date: 2023-12-05
#' date-modified: today
#' date-format: iso
#' format:
#'   html:
#'     toc: true
#'     toc-location: left
#'     toc-depth: 2
#'     number-sections: true
#'     smooth-scroll: true
#'     theme: readable
#'     code-copy: true
#'     code-download: true
#'     code-tools: true
#'     embed-resources: true
#'     anchor-sections: true
#'     html-math-method: katex
#' ---

#' # Introduction
#' 
#' This notebook contains several examples of how to use [Stan](https://mc-stan.org) in R with [__brms__](https://paul-buerkner.github.io/brms/). This notebook assumes basic knowledge of Bayesian inference and MCMC. The examples are related to [Bayesian data analysis course](https://avehtari.github.io/BDA_course_Aa/lto/).
#' 

#+ setup, include=FALSE
knitr::opts_chunk$set(cache=FALSE, message=FALSE, error=FALSE, warning=TRUE, comment=NA, out.width='95%')

#' **Load packages**
#| code-fold: true
#| cache: FALSE
library(tidyr)
library(dplyr)
library(tibble)
library(pillar)
library(stringr)
library(brms)
options(brms.backend = "cmdstanr", mc.cores = 2)
library(posterior)
options(posterior.num_args=list(digits=2))
options(pillar.negative = FALSE)
library(loo)
library(priorsense)
#options(priorsense.use_plot_theme=FALSE)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size=14))
library(tidybayes)
library(ggdist)
library(patchwork)
library(RColorBrewer)
library(tinytable)
options(tinytable_format_num_fmt = "significant_cell", tinytable_format_digits = 2, tinytable_tt_digits=2)
SEED <- 48927 # set random seed for reproducibility

#' # Bernoulli model
#' 
#' Toy data with sequence of failures (0) and successes (1). We would
#' like to learn about the unknown probability of success. brms wants
#' the data in a data frame format (tibble which is a variant of data
#' frame can be used, too).
data_bern <- data.frame(y = c(1, 1, 1, 0, 1, 1, 1, 0, 1, 0))

#' As usual in case of generalized linear models, (GLMs) brms defines
#' the priors on the latent model parameters. With Bernoulli the
#' default link function is logit, and thus the prior is set on
#' logit(theta). As there are no covariates logit(theta)=Intercept.
#' The brms default prior for Intercept is `student_t(3, 0, 2.5)`, but
#' we use `student_t(7, 0, 1.5)` which is close to logistic
#' distribution, and thus makes the prior near-uniform for theta.
#' We can simulate from these priors to check the implied prior on theta.
#' We next compare the result to using `normal(0, 1)` prior on logit
#' probability. We visualize the implied priors by sampling from the priors.
#'
#' `plogis` is the cumulative density function for logistic
#' distribution, which is equal to inverse-logit transformation.
#' Base-R's Student's t-distribution does not have `mu` and `sigma`
#' arguments, and thus we use a function from `ggdist` package.
data.frame(theta = plogis(ggdist::rstudent_t(n=20000, df=3, mu=0, sigma=2.5))) |>
  mcmc_hist() +
  xlim(c(0,1)) +
  labs(title='Default brms student_t(3, 0, 2.5) for Intercept')
data.frame(theta = plogis(ggdist::rstudent_t(n=20000, df=7, mu=0, sigma=1.5))) |>
  mcmc_hist() +
  xlim(c(0,1)) +
  labs(title='student_t(7, 0, 1.5) for Intercept')
#' Almost uniform prior on theta could be obtained also with normal(0,1.5)
data.frame(theta = plogis(rnorm(n=20000, mean=0, sd=1.5))) |>
  mcmc_hist() +
  xlim(c(0,1)) +
  labs(title='normal(0, 1.5) for Intercept')

#' Formula `y ~ 1` corresponds to a model $\mathrm{logit}(\theta) = \alpha\times 1 = \alpha$.
#' `brms` denotes the $\alpha$ as `Intercept`.
#| results: hide
fit_bern <- brm(y ~ 1, family = bernoulli(), data = data_bern,
                prior = prior(student_t(7, 0, 1.5), class='Intercept'),
                seed = SEED, refresh = 0)

#' Check the summary of the posterior and inference diagnostics.
fit_bern

#' Extract the posterior draws in data frame (df) format
draws <- as_draws_df(fit_bern)
#' We can select subset of stored variables with `subset_draws()`
#' and get summary information using `summarise_draws()`
draws |>
  subset_draws(variable='b_Intercept') |>
  summarise_draws()

#' We get a prettier table using `tinytable::tt()`. Earlier we had
#' defined options
#' ```
#' options(tinytable_format_num_fmt = "significant_cell", tinytable_format_digits = 2, tinytable_tt_digits=2)
#' ```
#' which makes the table to show only 2 significant digits for each value,
#' which is sufficient accuracy and having fewer digits makes the table
#' easier to read.
draws |>
  subset_draws(variable='b_Intercept') |>
  summarise_draws() |>
  tt()

#' We can compute the probability of success by using `plogis` which is
#' equal to inverse-logit function
draws <- draws |>
  mutate_variables(theta=plogis(b_Intercept))

#' Summary of `theta`
draws |>
  subset_draws(variable='theta') |>
  summarise_draws() |>
  tt()

#' Histogram of `theta` using `bayesplot::mcmc_hist()`
mcmc_hist(draws, pars='theta') +
  xlab('theta') +
  xlim(c(0,1))

#' Prior and likelihood sensitivity plot shows posterior density estimate
#' depending on amount of power-scaling. Overlapping lines indicate low
#' sensitivity and wider gaps between lines indicate greater sensitivity.
#| warning: false
powerscale_plot_dens(draws, fit=fit_bern, variable='theta',
                     help_text=FALSE) +
  xlim(c(0,1))

#' We can summarise the prior and likelihood sensitivity using
#' cumulative Jensen-Shannon distance.
powerscale_sensitivity(draws, fit=fit_bern, variable="theta") |>
  tt()

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
#| results: hide
fit_bin <- brm(y | trials(N) ~ 1, family = binomial(), data = data_bin,
               prior = prior(student_t(7, 0,1.5), class='Intercept'),
               seed = SEED, refresh = 0)

#' Check the summary of the posterior and inference diagnostics.
fit_bin

#' Extract the posterior draws in data frame format
draws <- as_draws_df(fit_bin)
#' Summary of latent intercept
draws |>
  subset_draws(variable='b_Intercept') |>
  summarise_draws() |>
  tt()

#' We can compute the probability of success by using `plogis()` which is
#' equal to inverse-logit function
draws <- draws |>
  mutate_variables(theta=plogis(b_Intercept))
#' Summary of theta
draws |>
  subset_draws(variable='theta') |>
  summarise_draws() |>
  tt(digits=2)

#' Histogram of `theta`
mcmc_hist(draws, pars='theta') +
  xlab('theta') +
  xlim(c(0,1))

#' Re-run the model with a new data dataset without recompiling using
#' argument `newdata`.
#| results: hide
#| cache: true
data_bin <- data.frame(N = c(5), y = c(4))
fit_bin <- update(fit_bin, newdata = data_bin)

#' Check the summary of the posterior and inference diagnostics.
fit_bin

#' Extract the posterior draws in data frame format
draws <- as_draws_df(fit_bin)
#' Summary of latent intercept
draws |>
  subset_draws(variable='b_Intercept') |>
  summarise_draws() |>
  tt(digits=2)

#' We can compute the probability of success by using `plogis()` which is
#' equal to inverse-logit function
draws <- draws |>
  mutate_variables(theta=plogis(b_Intercept))
#' Summary of `theta`
draws |>
  subset_draws(variable='theta') |>
  summarise_draws() |>
  tt(digits=2)

#' Histogram of `theta`
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
#' type (using `factor()` function), which is useful for categorical
#' variables.
data_bin2 <- data.frame(N = c(674, 680),
                        y = c(39,22),
                        grp2 = factor(c('control','treatment')))

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
#' treatment, and $x_\mathrm{treatment}$ is a vector with 1 for
#' treatment and 0 for control. As only of the vectors have 1, this
#' corresponds to separate models
#' $\mathrm{logit}(\theta_\mathrm{control}) = \beta_\mathrm{control}$
#' and $\mathrm{logit}(\theta_\mathrm{treatment}) =
#' \beta_\mathrm{treatment}$.  We can provide the same prior for all
#' $\beta$'s by setting the prior with `class='b'`. With prior
#' `student_t(7, 0,1.5)`, both $\beta$'s are shrunk towards 0, but
#' independently.
#| results: hide
fit_bin2 <- brm(y | trials(N) ~ 0 + grp2, family = binomial(),
                data = data_bin2,
                prior = prior(student_t(7, 0,1.5), class='b'),
                seed = SEED, refresh = 0)

#' Check the summary of the posterior and inference diagnostics. With `~ 0 +
#' grp2` there is no `Intercept` and \beta_\mathrm{control} and
#' \beta_\mathrm{treatment} are presented as `grp2control` and
#' `grp2treatment`.
fit_bin2

#' Compute theta for each group and the odds-ratio. `brms` uses
#' variable names `b_grp2control` and `b_grp2treatment` for
#' $\beta_\mathrm{control}$ and $\beta_\mathrm{treatment}$
#' respectively.
draws_bin2 <- as_draws_df(fit_bin2) |>
  mutate(theta_control = plogis(b_grp2control),
         theta_treatment = plogis(b_grp2treatment),
         oddsratio = (theta_treatment/(1-theta_treatment))/(theta_control/(1-theta_control)))
#' Plot histogram of odds-ratio
mcmc_hist(draws_bin2, pars='oddsratio') +
  scale_x_continuous(breaks=seq(0.2,1.6,by=0.2))+
  geom_vline(xintercept=1, linetype='dashed')

#' Compute probability that the oddsratio<1 and associated Monte Carlo
#' standard error (MCSE). To compute the probability we use comparison
#' `oddsratio<1` to find out for which posterior draws this is true,
#' and when we compute mean over TRUE and FALSE values, these values
#' are converted to 1 and 0, and the `mean()` is then equal to the
#' proportion of 1's. The usual MCSE estimate for mean can be use to
#' get the MCSE for this proportion.
draws_bin2 |>
  mutate(poddsratio = oddsratio<1) |>
  subset(variable='poddsratio') |>
  summarise_draws(mean, mcse_mean)  |>
  tt()

#' Compute odds-ratio 95% posterior interval, and associated MCSEs
draws_bin2 |>
  subset(variable='oddsratio') |>
  summarise_draws(~quantile(.x, probs = c(0.025, 0.975)), ~mcse_quantile(.x, probs = c(0.025, 0.975))) |>
  tt()


#' Make prior sensitivity analysis by power-scaling both prior and
#' likelihood.  Focus on odds-ratio which is the quantity of
#' interest. We see that the likelihood is much more informative than
#' the prior, and we would expect to see a different posterior only
#' with a highly informative prior (possibly based on previous similar
#' experiments).
#' Prior and likelihood sensitivity plot shows posterior density estimate
#' depending on amount of power-scaling. Overlapping lines indicate low
#' sensitivity and wider gaps between lines indicate greater sensitivity.
#| warning: false
powerscale_plot_dens(draws_bin2, fit=fit_bin2, variable='oddsratio',
                     help_text=FALSE) +
  labs(x='Odds-ratio', y=NULL) +
  scale_x_continuous(breaks=seq(0.2,1.4,by=0.2))+
  # reference line
  geom_vline(xintercept=1, linetype='dashed')

#' We can summarise the prior and likelihood sensitivity using
#' cumulative Jensen-Shannon distance. Here we prefer to show
#' 2 decimal digits (instead of the 2 significant digits we used before)
powerscale_sensitivity(draws_bin2, fit=fit_bin2, variable='oddsratio')  |>
  tt() |>
  format_tt(num_fmt="decimal")

#' Above we used formula `y | trials(N) ~ 0 + grp2` to have separate
#' model for control and treatment group. An alternative model `y |
#' trials(N) ~ grp2` which is equal to `y | trials(N) ~ 1 + grp2`,
#' would correspond to a model $\mathrm{logit}(\theta) = \alpha \times
#' 1 + \beta_\mathrm{treatment}\times x_\mathrm{treatment} = \alpha +
#' \beta_\mathrm{treatment}\times x_\mathrm{treatment}$. Now $\alpha$
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
  labs(title='student_t(7, 0, 1.5) for Intercept') +
  data.frame(theta_treatment = plogis(ggdist::rstudent_t(n=20000, df=7, mu=0, sigma=1.5))+
               plogis(ggdist::rstudent_t(n=20000, df=7, mu=0, sigma=1.5))) |>
  mcmc_hist() +
  xlim(c(0,1)) +
  labs(title='student_t(7, 0, 1.5) for Intercept and b_grp2treatment')

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
#' models the overall probability of death (via logistic link), and
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
#| results: hide
#| cache: true
fit_bin2 <- brm(y | trials(N) ~ 1 + (1 | grp2), family = binomial(),
                data = data_bin2,
                prior = prior(student_t(7, 0,1.5), class='Intercept'),
                seed = SEED, refresh = 0, control=list(adapt_delta=0.99))

#' Check the summary of the posterior and inference diagnostics. The summary
#' reports that there are Group-Level Effects: `~grp2` with 2 levels
#' (control and treatment), with `sd(Intercept)` denoting
#' $\sigma_\mathrm{grp}$. In addition, the summary lists
#' `Population-Level Effects: Intercept` ($\alpha$) as in the previous
#' non-hierarchical models.
fit_bin2

#' We can also look at the variable names `brms` uses internally.
#' We exclude variables `lprior` (log prior density) and `lp__`
#' (log posterior density).
as_draws_rvars(fit_bin2) |>
  subset_draws(variable=c('lprior','lp__'), exclude=TRUE) |>
  summarise_draws() |>
  tt()

#' Although there is no difference, illustrate how to compute the
#' oddsratio from a hierarchical model.
draws_bin2 <- as_draws_df(fit_bin2) |>
  mutate_variables(theta_control = plogis(b_Intercept + `r_grp2[control,Intercept]`),
                   theta_treatment = plogis(b_Intercept + `r_grp2[treatment,Intercept]`),
                   oddsratio = (theta_treatment/(1-theta_treatment))/(theta_control/(1-theta_control)))

draws_bin2 |> mcmc_hist(pars="oddsratio") +
  scale_x_continuous(breaks=seq(0.2,1.6,by=0.2))+
  geom_vline(xintercept=1, linetype='dashed')

#' Make also prior sensitivity analysis with focus on odds-ratio.
powerscale_sensitivity(draws_bin2, fit=fit_bin2, variable='oddsratio')  |>
  tt() |>
  format_tt(num_fmt="decimal")

#' # Linear Gaussian model
#' 
#' Use the Kilpisjärvi summer month temperatures 1952--2022 data from `aaltobda` package. We can read data directly from a URL.
load(url('https://github.com/avehtari/BDA_course_Aalto/raw/master/rpackage/data/kilpisjarvi2022.rda'))
data_lin <- data.frame(year = kilpisjarvi2022$year,
                       temp = kilpisjarvi2022$temp.summer)

#' Plot the data
data_lin |>
  ggplot(aes(year, temp)) +
  geom_point(color=2) +
  labs(x= "Year", y = 'Summer temp. @Kilpisjärvi') +
  guides(linetype = "none")

#' To analyse whether there has been change in the average summer month
#' temperature we use a linear model with Gaussian model for the
#' unexplained variation. By default brms uses uniform prior for the
#' coefficients.
#' 
#' Formula `temp ~ year` corresponds to model $\mathrm{temp} ~
#' \mathrm{normal}(\alpha + \beta \times \mathrm{temp}, \sigma)$.  The
#' model could also be defined as `temp ~ 1 + year` which explicitly
#' shows the intercept ($\alpha$) part. Using the variable names
#' `brms` uses the model can be written also as
#' $$
#' \mathrm{temp} \sim \mathrm{normal}(\mathrm{b\_Intercept}*1 + \mathrm{b\_year}*\mathrm{year}, \mathrm{sigma})
#' $$
#' We start with the default priors to see some tricks that `brms`
#' does behind the curtain.
#| results: hide
fit_lin <- brm(temp ~ year, data = data_lin, family = gaussian(),
               seed = SEED, refresh = 0)

#' Check the summary of the posterior and inference diagnostics.
fit_lin

#' Convergence diagnostics look good. We see that posterior mean of
#' `Intercept` is -34.7, which may sound strange, but that is the
#' intercept at year 0, that is, very far from the data range, and
#' thus doesn't have meaningful interpretation directly. The posterior
#' mean of `year` coefficient is 0.02, that is, we estimate that the
#' summer temperature is increasing 0.02°C per year (which would make
#' 1°C in 50 years).
#' 
#' We can check $R^2$ which corresponds to the proportion of variance
#' explained by the model. The linear model explains 0.16=16% of the
#' total data variance.
bayes_R2(fit_lin) |>
  as_tibble() |>
  tt()

#' We can check the all the priors used with `prior_summary()`
prior_summary(fit_lin) |>
  tt()

#' We see that `class=b` and `coef=year` have prior `flat`, that is,
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
#' covariate values.
#'
#' We can turn of the automatic centering by replacing the formula
#' with `bf(temp ~ year, center=FALSE)`. Or we can set the prior on
#' original intercept by using a formula `temp ~ 0 + Intercept +
#' year`.
#'
#' In this case, we are
#' happy with the default prior for the intercept. In this specific
#' case, the flat prior on coefficient is also fine, but we add an
#' weakly informative prior just for the illustration. Let's assume we
#' expect the temperature to change less than 1°C in 10 years. With
#' `student_t(3, 0, 0.03)` about 95% prior mass has less than 0.1°C
#' change in year, and with low degrees of freedom (3) we have thick
#' tails making the likelihood dominate in case of prior-data
#' conflict. In real life, we do have much more information about the
#' temperature change, and naturally a hierarchical spatio-temporal
#' model with all temperature measurement locations would be even
#' better.
#| results: hide
fit_lin <- brm(temp ~ year, data = data_lin, family = gaussian(),
               prior = prior(student_t(3, 0, 0.03), class='b'),
               seed = SEED, refresh = 0)

#' Check the summary of the posterior and inference diagnostics.
fit_lin

#' Make prior sensitivity analysis by power-scaling both prior and likelihood.
powerscale_sensitivity(fit_lin)  |>
  tt() |>
  format_tt(num_fmt="decimal")

#' Our weakly informative proper prior has negligible sensitivity, and
#' the likelihood is informative.

#' Extract the posterior draws and check the summaries. We exclude
#' variables `lprior` (log prior density) and `lp__` (log posterior
#' density).
draws_lin <- as_draws_df(fit_lin) 
draws_lin |>
  subset_draws(variable=c('lprior','lp__'), exclude=TRUE) |>
  summarise_draws() |> 
  tt()

#' Histogram of `b_year`
draws_lin |>
  mcmc_hist(pars='b_year') +
  xlab('Average temperature increase per year')

#' Compute the probability that the coefficient `b_year` > 0 and the
#' corresponding MCSE.
draws_lin |>
  mutate(I_b_year_gt_0 = b_year>0) |>
  subset_draws(variable='I_b_year_gt_0') |>
  summarise_draws(mean, mcse_mean)  |>
  tt(digits=2)
#' All posterior draws have `b_year>0`, the probability gets rounded
#' to 1, and MCSE is not available as the observed posterior variance
#' is 0.
#' 

#' 95% posterior interval for temperature increase per 100 years and
#' associated MCSEs.
draws_lin |>
  mutate(b_year_100 = b_year*100) |>
  subset_draws(variable='b_year_100') |>
  summarise_draws(~quantile(.x, probs = c(0.025, 0.975)),
                  ~mcse_quantile(.x, probs = c(0.025, 0.975))) |>
  tt() |>
  format_tt(num_fmt="decimal")

#' Plot posterior draws of the linear function values at each year.
#' `add_linpred_draws()` takes the years from the data passed via
#' pipe, and uses `fit_lin` to make the linear model predictions.
#' `add_linpred_draws()` corresponds to `brms::posterior_linpred()`
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

#' Plot a spaghetti plot for 100 posterior draws.
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

#' Plot posterior predictive distribution at each year until 2030.
#' `add_predicted_draws()` takes the years from the data and uses
#' `fit_lin` to draw from the posterior predictive distribution.
#' `add_predicted_draws()` corresponds to `brms::posterior_predict()`.
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

#' Posterior predictive check with density overlays examines the whole
#' temperature distribution. We generate replicate data using 20 different
#' posterior draws (with argument `ndraws`).
pp_check(fit_lin, type='dens_overlay', ndraws=20)

#' LOO-PIT check is good for checking whether the normal distribution
#' is well describing the variation as it examines the calibration
#' of LOO predictive distributions conditionally on each year. LOO-PIT
#' plot looks good. We use all posterior draws to estimate LOO predictive
#' distributions.
pp_check(fit_lin, type='loo_pit_qq', ndraws=4000)

#' # Linear Student's $t$ model
#' 
#' The temperatures used in the above analyses are averages over three
#' months, which makes it more likely that they are normally
#' distributed, but there can be extreme events in the feather and we
#' can check whether more robust Student's $t$ observation model would
#' give different results (although LOO-PIT check did already indicate
#' that the normal would be good).
#' 
#+  results='hide'
fit_lin_t <- brm(temp ~ year, data = data_lin,
                 family = student(),
                 prior = prior(student_t(3, 0, 0.03), class='b'),
                 seed = SEED, refresh = 0)

#' Check the summary of the posterior and inference diagnostics. The
#' `b_year` posterior looks similar as before and the posterior for
#' degrees of freedom `nu` has most of the posterior mass for quite
#' large values indicating there is no strong support for thick tailed
#' variation in average summer temperatures.
fit_lin_t

#' # Pareto-smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO)
#' 
#' We can use leave-one-out cross-validation to compare the expected
#' predictive performance.
#' 
#' LOO comparison shows normal and Student's $t$ model have similar
#' performance. As `loo_compare()` returns it's own specific object
#' type, we need to do some manipulation to change it to a data frame
#' suitable for `tt()`.
loo_compare(loo(fit_lin), loo(fit_lin_t)) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  select(model, elpd_diff, se_diff) |>
  tt()

#' # Heteroskedastic linear model
#' 
#' Heteroscedasticity assumes that the variation around the linear
#' mean can also vary. We can allow sigma to depend on year, too.
#' `brms` supports multiple models using `bf()` function, and `sigma`
#' is a special keyword for the residual standard deviation.
#' Although the additional component is written as `sigma ~ year`, the
#' log link function is used and the model is for log(sigma). 
#+  results='hide'
fit_lin_h <- brm(bf(temp ~ year,
                    sigma ~ year),
                 data = data_lin, family = gaussian(),
                 prior = prior(student_t(3, 0, 0.03), class='b'),
                 seed = SEED, refresh = 0)

#' Check the summary of the posterior and inference diagnostics. The `b_year`
#' posterior looks similar as before. The posterior for `sigma_year`
#' looks like having most of the mass for negative values, indicating
#' decrease in temperature variation around the mean.
fit_lin_h

#' Histograms of `b_year` and `b_sigma_year`
as_draws_df(fit_lin_h) |>
  mcmc_areas(pars=c('b_year', 'b_sigma_year'))

#' As $log(x)$ is almost linear when $x$ is close to zero, we can see
#' that the `sigma` is decreasing about 1% per year (95% interval from
#' 0% to 2%).
#' 

#' Plot the posterior predictive distribution at each year until 2030
#' `add_predicted_draws()` takes the years from the data and uses
#' `fit_lin_h` to draw from the posterior predictive distribution.
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

#' Make prior sensitivity analysis by power-scaling both prior and likelihood.
powerscale_sensitivity(fit_lin_h)  |>
  tt() |>
  format_tt(num_fmt="decimal")

#' We can use leave-one-out cross-validation to compare the expected
#' predictive performance.
#' 
#' LOO comparison shows homoskedastic normal and heteroskedastic
#' normal models have similar performances.
loo_compare(loo(fit_lin), loo(fit_lin_h)) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  select(model, elpd_diff, se_diff) |>
  tt()

#' # Heteroskedastic non-linear model
#' 
#' We can test the linearity assumption by using non-linear spline
#' functions, by using `s(year)` terms. Sampling is slower as the
#' posterior gets more complex.
#' 
#+  results='hide'
fit_spline_h <- brm(bf(temp ~ s(year),
                       sigma ~ s(year)),
                    data = data_lin, family = gaussian(),
                    seed = SEED, refresh = 0)

#' We get warnings about divergences, and try rerunning with higher
#' adapt_delta, which leads to using smaller step sizes. Often
#' `adapt_delta=0.999` leads to very slow sampling and is not
#' generally recommended, but with this small data, this is not an
#' issue.
#| results: hide
#| cache: true
fit_spline_h <- update(fit_spline_h, control = list(adapt_delta=0.999))

#' Check the summary of the posterior and inference diagnostics. We're not
#' anymore able to make interpretation of the temperature increase
#' based on this summary. For splines, we see prior scales `sds` for
#' the spline coefficients.
fit_spline_h

#' We can still plot posterior predictive distribution at each year
#' until 2030 `add_predicted_draws()` takes the years from the data
#' and uses `fit_lin_h` to draw from the posterior predictive distribution.
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

#' And we can use LOO-CV to compare the
#' expected predictive performance.
#' LOO comparison shows homoskedastic normal linear and
#' heteroskedastic normal spline models have similar
#' performances. There are not enough observations to make clear
#' difference between the models.
loo_compare(loo(fit_lin), loo(fit_spline_h)) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  select(model, elpd_diff, se_diff) |>
  tt()

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
                  ~mcse_quantile(.x, probs = c(0.025, 0.975))) |>
  tt() |>
  format_tt(num_fmt="decimal")

#' Make prior sensitivity analysis by power-scaling both prior and
#' likelihood with focus on average summer temperature increase from
#' 1952 to 2022.
powerscale_sensitivity(temp_diff, fit=fit_spline_h, variable='temp_diff') |>
  tt() |>
  format_tt(num_fmt="decimal")

#' Probability that the average summer temperature has increased from
#' 1952 to 2022 is 99.5%.
temp_diff |>
  mutate(I_temp_diff_gt_0 = temp_diff>0,
         temp_diff = NULL) |>
  subset_draws(variable='I_temp_diff_gt_0') |>
  summarise_draws(mean, mcse_mean) |>
  tt(digits=2)


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
#| results: hide
#| cache: true
fit_pooled <- brm(quality ~ 1, data = factory, refresh=0)

#' Check the summary of the posterior and inference diagnostics.
fit_pooled

#' ## Separate model
#'
#' As comparison make also separate model. To make it completely
#' separate we need to have different sigma for each machine, too.
#| results: hide
#| cache: true
fit_separate <- brm(bf(quality ~ 0 + machine,
                       sigma ~ 0 + machine),
                    data = factory, refresh=0)

#' Check the summary of the posterior and inference diagnostics.
fit_separate

#' # Common variance hierarchical model (ANOVA)
#| results: hide
#| cache: true
fit_hier <- brm(quality ~ 1 + (1 | machine),
                data = factory, refresh = 0)

#' Check the summary of the posterior and inference diagnostics.
fit_hier

#' LOO comparison shows the hierarchical model is the best. The
#' differences are small as the number of observations is small and
#' there is a considerable prediction (aleatoric) uncertainty.
loo_compare(loo(fit_pooled), loo(fit_separate), loo(fit_hier)) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  select(model, elpd_diff, se_diff) |>
  tt()

#' Different model posterior distributions for the mean
#' quality. Pooled model ignores the variation between
#' machines. Separate model doesn't take benefit from the similarity of
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

#' Make prior sensitivity analysis by power-scaling both prior and
#' likelihood with focus on mean quality of each machine. We see no
#' prior sensitivity.
machine_mean <- fit_hier |>
  as_draws_df() |>
  mutate(across(matches('r_machine'), ~ .x - b_Intercept)) |>
  subset_draws(variable='r_machine', regex=TRUE) |>
  set_variables(paste0('machine_mean[', 1:6, ']'))
powerscale_sensitivity(machine_mean, fit=fit_hier, variable='machine_mean') |>
  tt() |>
  format_tt(num_fmt="decimal")

#' 
#' # Hierarchical binomial model
#'
#' [Sorafenib Toxicity Dataset in `metadat` R package](https://wviechtb.github.io/metadat/reference/dat.ursino2021.html)
#' includes results from 13 studies investigating the occurrence of
#' dose limiting toxicities (DLTs) at different doses of Sorafenib.
#'
#' Load data
load(url('https://github.com/wviechtb/metadat/raw/master/data/dat.ursino2021.rda'))
head(dat.ursino2021)

#' Number of patients per study
dat.ursino2021 |>
  group_by(study) |>
  summarise(N = sum(total)) |>
  ggplot(aes(x=N, y=study)) +
  geom_col(fill=4) +
  labs(x='Number of patients per study', y='Study')

#' Distribution of doses
dat.ursino2021 |>
  ggplot(aes(x=dose)) +
  geom_histogram(breaks=seq(50,1050,by=100), fill=4, colour=1) +
  labs(x='Dose (mg)', y='Count') +
  scale_x_continuous(breaks=seq(100,1000,by=100))

#' Each study is using 2--6 different dose levels. Three studies
#' that include only two dose levels are likely to provide weak
#' information on slope.
crosstab <- with(dat.ursino2021,table(dose,study))
data.frame(count=colSums(crosstab), study=colnames(crosstab)) |>
  ggplot(aes(x=count, y=study)) +
  geom_col(fill=4) +
  labs(x='Number of dose levels per study', y='Study')

#' Pooled model assumes all studies have the same dose effect
#' (reminder: `~ dose` is equivalent to `~ 1 + dose`).
#' We use similar priors as in earlier binomial models.
#| results: hide
#| cache: true
fit_pooled <- brm(events | trials(total) ~ dose,
                  prior = c(prior(student_t(7, 0, 1.5), class='Intercept'),
                            prior(normal(0, 1), class='b')),
                  family=binomial(), data=dat.ursino2021)

#' Check the summary of the posterior and inference diagnostics.
fit_pooled

#' Dose coefficient seems to be very small. Looking at the posterior,
#' we see that it is positive with high probability. 
fit_pooled |>
  as_draws() |>
  subset_draws(variable='b_dose') |>
  summarise_draws(~quantile(.x, probs = c(0.025, 0.975)), ~mcse_quantile(.x, probs = c(0.025, 0.975))) |>
  tt()

#' The dose was reported in mg, and most values are in hundreds. It is
#' often sensible to switch to a scale in which the range of values is
#' closer to unit range. In this case it is natural to use g instead
#' of mg.
dat.ursino2021 <- dat.ursino2021 |>
  mutate(doseg = dose/1000)

#' Fit the pooled model again using `doseg`
#| results: hide
#| cache: true
fit_pooled <- brm(events | trials(total) ~ doseg,
                  prior = c(prior(student_t(7, 0, 1.5), class='Intercept'),
                            prior(normal(0, 1), class='b')),
                  family=binomial(), data=dat.ursino2021)

#' Check the summary of the posterior and inference diagnostics.
fit_pooled

#' Now it is easier to interpret the presented values.
#'
#' Prior and likelihood sensitivity plot shows posterior density estimate
#' depending on amount of power-scaling. Overlapping lines indicate low
#' sensitivity and wider gaps between lines indicate greater sensitivity.
#| warning: false
fit_pooled |>
  powerscale_plot_dens(variable='b_doseg', help_text=FALSE) +
  labs(x='Dose (g) coefficient', y=NULL) 

#' We see a strong prior prior sensitivity.
#'
#' Power-scaling with cumulative Jensen-Shannon distance diagnostic
#' indicates prior-data conflict.
powerscale_sensitivity(fit_pooled, variable='b_doseg') |>
  tt() |>
  format_tt(num_fmt="decimal")

#' Comparing the posterior of `b_doseg` (90\%-interval [1.3, 3.6]) to
#' the prior normal(0,1), we see that when we scaled the covariate, we
#' forgot to check that the prior still makes sense.
#'
#' We make prior predictive checking by fixing the intercept=0 (or
#' otherwise prior variation from the intercept would hide the prior
#' variation from the `b_doseg` and centering the covariate (this is
#' what `brms` does by default).
data.frame(theta=plogis(sample(with(dat.ursino2021, doseg - mean(doseg)),
                               4000,replace=TRUE) *
                          rnorm(4000, mean=0, sd=1))) |>
  ggplot(aes(x=theta)) +
  stat_slab() +
  xlim(c(0,1)) +
  scale_y_continuous(breaks=NULL) +
  ylab('') +
  theme(axis.line.y=element_blank())

#' We see that the prior is very informative (here the width of the
#' prior matters, as the location would be based on the intercept.
#'
#' We could have made the prior predictive checking earlier, but here
#' we intentionally wanted to illustrate how `priorsense` can catch
#' prior problems.
#' 
#' Checking that sd of `doseg` is about 1/5, 
sd(dat.ursino2021$doseg) |> round(2)
#' We adjust the prior to be normal(0, 5), so that expected sd of
#' `doseg * b_doseg` is about 1. Prior predictive checking looks better now.
data.frame(theta=plogis(sample(with(dat.ursino2021, doseg - mean(doseg)),
                               4000,replace=TRUE) *
                          rnorm(4000, mean=0, sd=5))) |>
  ggplot(aes(x=theta)) +
  stat_slab() +
  xlim(c(0,1)) +
  scale_y_continuous(breaks=NULL) +
  ylab('') +
  theme(axis.line.y=element_blank())

#' The narrower part in the prior distribution is due to the data
#' having more small doses than large doses.
#'
#' We refit the model with the wider prior.
#| results: hide
#| cache: true
fit_pooled <- brm(events | trials(total) ~ doseg,
                  prior = c(prior(student_t(7, 0, 1.5), class='Intercept'),
                            prior(normal(0, 5), class='b')),
                  family=binomial(), data=dat.ursino2021)
fitp_pooled <- update(fit_pooled, sample_prior='only')

#' And the prior-data conflict has gone.
#| warning: false
fit_pooled |>
  powerscale_plot_dens(variable='b_doseg', help_text=FALSE) +
  labs(x='Dose (g) coefficient', y=NULL)

powerscale_sensitivity(fit_pooled, variable='b_doseg') |>
  tt() |>
  format_tt(num_fmt="decimal")

#' Separate model assumes all studies have different dose effect.
#' It would be a bit complicated to set a different prior on study specific
#' intercepts and other coefficients, so we use the same prior for all.
#| results: hide
#| cache: true
fit_separate <- brm(events | trials(total) ~ 0 + study + doseg:study,
                    prior=prior(student_t(7, 0, 5), class='b'),
                    family=binomial(), data=dat.ursino2021)

#' Check the summary of the posterior and inference diagnostics.
fit_separate

#' We build two different hierarchical models. The first one has
#' hierarchical model for the intercept, that is, each study has a
#' parameter telling how much that study differs from the common
#' population intercept.
#| results: hide
#| cache: true
fit_hier1 <- brm(events | trials(total) ~ doseg + (1 | study),
                 family=binomial(),
                 prior=c(prior(student_t(7, 0, 3), class='Intercept'),
                         prior(normal(0, 5), class='b')),
                 data=dat.ursino2021)

#' The second hierarchical model assumes that also the slope can vary
#' between the studies.
#| results: hide
#| cache: true
fit_hier2 <- brm(events | trials(total) ~ doseg + (doseg | study),
                 family=binomial(),
                 prior=c(prior(student_t(7, 0, 10), class='Intercept'),
                         prior(normal(0, 5), class='b')),
                 data=dat.ursino2021)

#' We seem some divergences due to highly varying posterior
#' curvature. We repeat the sampling with higher adapt_delta, which
#' adjust the step size to be smaller. Higher adapt_delta makes the
#' computation slower, but that is not an issue in this case. If you
#' get divergences with `adapt_delta=0.99`, it is likely that even
#' larger values don't help, and you need to consider different
#' parameterization, different model, or more informative priors.
#| results: hide
#| cache: true
fit_hier2 <- update(fit_hier2, control=list(adapt_delta=0.99))

#' LOO-CV comparison
loo_compare(loo(fit_pooled), loo(fit_separate), loo(fit_hier1), loo(fit_hier2)) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  select(model, elpd_diff, se_diff) |>
  tt()

#' We get warnings about several Pareto k's > 0.7 in PSIS-LOO for
#' separate model, but as in that case the LOO-CV estimate is usually
#' overoptimistic and the separate model is the worst, there is no
#' need to use more accurate computation for the separate model.
#'
#' We get warnings about a few Pareto k's > 0.7 in PSIS-LOO for both
#' hierarchical models. We can improve the accuracy be running MCMC
#' for these LOO folds. We use `add_criterion()` function to store the
#' LOO computation results as they take a bit longer now. We get some
#' divergences in case of the second hierarchical model, as leaving
#' out an observation for a study that has only two dose levels is
#' making the posterior having a difficult shape.
#| results: hide
#| cache: true
fit_hier1 <- add_criterion(fit_hier1, criterion='loo', reloo=TRUE)
fit_hier2 <- add_criterion(fit_hier2, criterion='loo', reloo=TRUE)

#' We repeat the LOO-CV comparison (without separate model). `loo()`
#' function is using the results added to the fit objects.
loo_compare(loo(fit_pooled), loo(fit_hier1), loo(fit_hier2)) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  select(model, elpd_diff, se_diff) |>
  tt()

#' The results did not change much. The first hierarchical model is
#' slightly better than other models, but for predictive purposes
#' there is not much difference (there is high aleatoric uncertainty
#' in the predictions). Adding hierarchical model for the slope,
#' decreased the predictive performance and thus it is likely that
#' there is not enough information about the variation in slopes
#' between studies.
#' 

#' Posterior predictive checking showing the observed and predicted
#' number of events. Rootogram uses square root of counts on y-axis for
#' better scaling. Rootogram is useful for count data when the range
#' of counts is small or moderate.
pp_check(fit_pooled, type = "rootogram") +
  labs(title='Pooled model')
pp_check(fit_hier1, type = "rootogram") +
  labs(title='Hierarchical model 1')
pp_check(fit_hier2, type = "rootogram") +
  labs(title='Hierarchical model 2')

#' We see that the hierarchical models have higher probability for
#' future counts that are bigger than maximum observed count and
#' longer predictive distribution tail. This is natural as uncertainty
#' in the variation between studies increases predictive uncertainty,
#' too, especially as the number of studies is relatively small.
#' 

#' The population level coefficient posterior given pooled model
plot_posterior_pooled <- mcmc_areas(as_draws_df(fit_pooled), regex_pars='b_doseg') +
  geom_vline(xintercept=0, linetype='dashed') +
  labs(title='Pooled model')
#' The population level coefficient posterior given hierarchical model 1
plot_posterior_hier1 <- mcmc_areas(as_draws_df(fit_hier1), regex_pars='b_doseg') +
  geom_vline(xintercept=0, linetype='dashed') +
  labs(title='Hierarchical model 1')
#' The population level coefficient posterior given hierarchical model 3
plot_posterior_hier2 <- mcmc_areas(as_draws_df(fit_hier2), regex_pars='b_doseg') +
  geom_vline(xintercept=0, linetype='dashed') +
  labs(title='Hierarchical model 2')

(plot_posterior_pooled / plot_posterior_hier1 / plot_posterior_hier2) * xlim(c(0,8.5))

#' All models agree that the slope is very likely positive. The
#' hierarchical models have more uncertainty, but also higher
#' posterior mean.
#' 

#' When we look at the study specific parameters, we see that the
#' Miller study has slightly higher intercept (leading to higher theta).
(mcmc_areas(as_draws_df(fit_hier1), regex_pars='r_study\\[.*Intercept') +
   labs(title='Hierarchical model 1')) /
  (mcmc_areas(as_draws_df(fit_hier2), regex_pars='r_study\\[.*Intercept') +
     labs(title='Hierarchical model 2'))

#'
#' There are no clear differences in slopes.
mcmc_areas(as_draws_df(fit_hier2), regex_pars='r_study\\[.*doseg') +
  labs(title='Hierarchical model 2')
#'

#' Based on LOO comparison we could continue with any of the models
#' except the separate one, but if we want to take into account the
#' unknown possible study variations, it is best to continue with the
#' hierarchical model 2.  We could reduce the uncertainty by spending
#' some effort to elicit a more informative priors for the between
#' study variation, by searching open study databases for similar
#' studies. In this example, we skip that and continue with other
#' parts of the workflow.

#' We check the prior sensitivity in hierarchical model 2
#| warning: false
fit_hier2 |>
  powerscale_plot_dens(variable='b_doseg', help_text=FALSE) +
  labs(x='Dose (g) coefficient', y=NULL)

powerscale_sensitivity(fit_hier2, variable='b_doseg') |>
  tt() |>
  format_tt(num_fmt="decimal")

#' 
#' The posterior for the probability of event given certain dose and a
#' new study for hierarchical model 2.
data.frame(study='new',
           doseg=seq(0.1,1,by=0.1),
           dose=seq(100,1000,by=100),
           total=1) |>
  add_linpred_draws(fit_hier2, transform=TRUE, allow_new_levels=TRUE) |>
  ggplot(aes(x=dose, y=.linpred)) +
  stat_lineribbon(.width = c(.95), alpha = 1/2, color=brewer.pal(5, "Blues")[[5]]) +
  scale_fill_brewer()+
  labs(x= "Dose (mg)", y = 'Probability of event', title='Hierarchical model') +
  theme(legend.position="none") +
  geom_hline(yintercept=0) +
  scale_x_continuous(breaks=seq(100,1000,by=100)) 

#' If we plot individual posterior draws, we see that there is a lot of
#' uncertainty about the overall probability (explained by the
#' variation in Intercept in different studies), but less uncertainty
#' about the slope.
data.frame(study='new',
           doseg=seq(0.1,1,by=0.1),
           dose=seq(100,1000,by=100),
           total=1) |>
  add_linpred_draws(fit_hier2, transform=TRUE, allow_new_levels=TRUE, ndraws=100) |>
  ggplot(aes(x=dose, y=.linpred)) +
  geom_line(aes(group=.draw), alpha = 1/2, color = brewer.pal(5, "Blues")[[3]])+
  scale_fill_brewer()+
  labs(x= "Dose (g)", y = 'Probability of event') +
  theme(legend.position="none") +
  geom_hline(yintercept=0) +
  scale_x_continuous(breaks=seq(100,1000,by=100)) 

#' 
#' # Hierarchical binomial model 2
#'
#' [Studies on Pharmacologic Treatments for Chronic Obstructive
#' Pulmonary
#' Disease](https://wviechtb.github.io/metadat/reference/dat.baker2009.html)
#' includes results from 39 trials examining pharmacologic treatments
#' for chronic obstructive pulmonary disease (COPD).
#'
#' Load data
load(url('https://github.com/wviechtb/metadat/raw/master/data/dat.baker2009.rda'))
# force character strings to factors for easier plotting
dat.baker2009 <- dat.baker2009 |>
  mutate(study = factor(study),
         treatment = factor(treatment),
         id = factor(id))

#' Look at six first lines of the data frame
head(dat.baker2009)

#' Total number of patients in each study varies a lot
dat.baker2009 |>
  group_by(study) |>
  summarise(N = sum(total)) |>
  ggplot(aes(x=N, y=study)) +
  geom_col(fill=4) +
  labs(x='Number of patients per study', y='Study')

#' None of the treatments is included in every study, and each study includes
#' $2--4$ treatments.
crosstab <- with(dat.baker2009,table(study, treatment))
#
plot_treatments <- data.frame(number_of_studies=colSums(crosstab), treatment=colnames(crosstab)) |>
  ggplot(aes(x=number_of_studies,y=treatment)) +
  geom_col(fill=4) +
  labs(x='Number of studies with a treatment X', y='Treatment') +
  geom_vline(xintercept=nrow(crosstab), linetype='dashed') +
  scale_x_continuous(breaks=c(0,10,20,30,39))
#
plot_studies <- data.frame(number_of_treatments=rowSums(crosstab), study=rownames(crosstab)) |>
  ggplot(aes(x=number_of_treatments,y=study)) +
  geom_col(fill=4) +
  labs(x='Number of treatments in a study Y', y='Study') +
  geom_vline(xintercept=ncol(crosstab), linetype='dashed') +
  scale_x_continuous(breaks=c(0,2,4,6,8))
#
plot_treatments + plot_studies

#' The following plot shows which treatments were in which studies.
crosstab |>
  as_tibble() |>
  ggplot(aes(x=study, y=treatment, fill=as.factor(n))) +
  geom_tile() +
  scale_fill_manual(name = "", values = c("white",4)) +
  labs(x="Study", y="Treatment") +
  theme(aspect.ratio=2*8/39,
        legend.position='none',
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

#' The first model is pooling the information over studies, but estimating separate
#' theta for each treatment (including placebo).
#| results: hide
#| cache: true
fit_pooled <- brm(exac | trials(total) ~ 0 + treatment,
                  prior = prior(student_t(7, 0, 1.5), class='b'),
                  family=binomial(), data=dat.baker2009)

#' Check the summary of the posterior and inference diagnostics.
fit_pooled

#' Treatment effect posteriors
fit_pooled |>
  as_draws_df() |>
  subset_draws(variable='b_', regex=TRUE) |>
  set_variables(paste0('b_treatment[', levels(factor(dat.baker2009$treatment)), ']')) |>
  as_draws_rvars() |>
  spread_rvars(b_treatment[treatment]) |>
  mutate(theta_treatment = rfun(plogis)(b_treatment)) |>
  ggplot(aes(xdist=theta_treatment, y=treatment)) +
  stat_halfeye() +
  labs(x='theta', y='Treatment', title='Pooled over studies, separate over treatments')  

#' Treatment effect odds-ratio posteriors
theta <- fit_pooled |>
  as_draws_df() |>
  subset_draws(variable='b_', regex=TRUE) |>
  set_variables(paste0('b_treatment[', levels(factor(dat.baker2009$treatment)), ']')) |>
  as_draws_rvars() |>
  spread_rvars(b_treatment[treatment]) |>
  mutate(theta_treatment = rfun(plogis)(b_treatment))
theta_placebo <- filter(theta,treatment=='Placebo')$theta_treatment[[1]]
theta |>
  mutate(treatment_oddsratio = (theta_treatment/(1-theta_treatment))/(theta_placebo/(1-theta_placebo))) |>
  filter(treatment != "Placebo") |>
  ggplot(aes(xdist=treatment_oddsratio, y=treatment)) +
  stat_halfeye() +
  labs(x='Odds-ratio', y='Treatment', title='Pooled over studies, separate over treatments') +
  geom_vline(xintercept=1, linetype='dashed')

#' We see a big variation between treatments and two treatments seem
#' to be harmful, which is suspicious. Looking at the data we see that
#' not all studies included all treatments, and thus if some of the
#' studies had more events, then the above estimates can be wrong.
#' 

#' The target is discrete count, but as the range of counts is big, a
#' rootogram would look messy, and density overlay plot is a better
#' choice.  Posterior predictive checking with kernel density
#' estimates for the data and 10 posterior predictive replicates shows
#' clear discrepancy.
pp_check(fit_pooled, type='dens_overlay')

#' Posterior predictive checking with PIT values and ECDF difference
#' plot with envelope shows clear discrepancy.
pp_check(fit_pooled, type='pit_ecdf', ndraws=4000)

#' Posterior predictive checking with LOO-PIT values show clear discrepancy.
pp_check(fit_pooled, type='loo_pit_qq', ndraws=4000) +
  geom_abline() +
  ylim(c(0,1))

#' The second model uses a hierarchical model both for treatment
#' effects and study effects.
#| results: hide
#| cache: true
fit_hier <- brm(exac | trials(total) ~ (1 | treatment) + (1 | study),
                family=binomial(), data=dat.baker2009)

#' Check the summary of the posterior and inference diagnostics.
fit_hier

#' LOO-CV comparison
loo_compare(loo(fit_pooled), loo(fit_hier)) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  select(model, elpd_diff, se_diff) |>
  tt()

#' We get warnings about Pareto k's > 0.7 in PSIS-LOO, but as the
#' difference between the models is huge, we can be confident that the
#' order would the same if we fixed the computation, and the
#' hierarchical model is much better and there is high variation
#' between studies. Clearly there are many highly influential
#' observations.
#'
#' Posterior predictive checking with kernel density estimates for the
#' data and 10 posterior predictive replicates looks good (although
#' with this many parameters, this check is likely to be optimistic).
pp_check(fit_hier, type='dens_overlay')

#' Posterior predictive checking with PIT values and ECDF difference
#' plot with envelope looks good (although with this many parameters,
#' this check is likely to be optimistic).
pp_check(fit_hier, type='pit_ecdf', ndraws=4000)

#' Posterior predictive checking with LOO-PIT values look good
#' (although as there are Pareto-khat warnings, it is possible that
#' this diagnostic is optimistic).
pp_check(fit_hier, type='loo_pit_qq', ndraws=4000) +
  geom_abline() +
  ylim(c(0,1))

#' Treatment effect posteriors have now much less variation.
fit_hier |>
  spread_rvars(b_Intercept, r_treatment[treatment,]) |>
  mutate(theta_treatment = rfun(plogis)(b_Intercept + r_treatment)) |>
  ggplot(aes(xdist=theta_treatment, y=treatment)) +
  stat_halfeye() +
  labs(x='theta', y='Treatment', title='Hierarchical over studies, hierarchical over treatments')  

#' Study effect posteriors show the expected high variation.
fit_hier |>
  spread_rvars(b_Intercept, r_study[study,]) |>
  mutate(theta_study = rfun(plogis)(b_Intercept + r_study)) |>
  ggplot(aes(xdist=theta_study, y=study)) +
  stat_halfeye() +
  labs(x='theta', y='Study', title='Hierarchical over studies, hierarchical over treatments')  

#' Treatment effect odds-ratio posteriors
theta <- fit_hier |>
  spread_rvars(b_Intercept, r_treatment[treatment,]) |>
  mutate(theta_treatment = rfun(plogis)(b_Intercept + r_treatment))
theta_placebo <- filter(theta,treatment=='Placebo')$theta_treatment[[1]]
theta |>
  mutate(treatment_oddsratio = (theta_treatment/(1-theta_treatment))/(theta_placebo/(1-theta_placebo))) |>
  filter(treatment != "Placebo") |>
  ggplot(aes(xdist=treatment_oddsratio, y=treatment)) +
  stat_halfeye() +
  labs(x='Odds-ratio', y='Treatment', title='Hierarchical over studies, hierarchical over treatments') +
  geom_vline(xintercept=1, linetype='dashed')

#' Treatment effect odds-ratios look now more reasonable. As now all
#' treatments were compared to placebo, there is less overlap in the
#' distributions as when looking at the thetas, as all thetas include
#' similar uncertainty about the overall theta due to high variation
#' between studies.

#' We check the prior sensitivity with focus in odds-ratios in hierarchical model
theta <- fit_hier |>
  spread_rvars(b_Intercept, r_treatment[treatment,]) |>
  mutate(theta_treatment = rfun(plogis)(b_Intercept + r_treatment))
theta_placebo <- filter(theta,treatment=='Placebo')$theta_treatment[[1]]
oddsratio <- theta |>
  mutate(treatment_oddsratio = (theta_treatment/(1-theta_treatment))/(theta_placebo/(1-theta_placebo))) |>
  filter(treatment != "Placebo")
oddsratio <- oddsratio$treatment_oddsratio |>
  as_draws_df() |>
  set_variables(oddsratio$treatment)
powerscale_sensitivity(oddsratio, fit=fit_hier) |>
  tt() |>
  format_tt(num_fmt="decimal")

#' The third model includes interaction so that the treatment can depend on study. 
#| results: hide
#| cache: true
fit_hier2 <- brm(exac | trials(total) ~ (1 | treatment) + (treatment | study),
                 family=binomial(), data=dat.baker2009, control=list(adapt_delta=0.9))

#' LOO comparison shows that the added interaction doesn't improve the model.
loo_compare(loo(fit_hier), loo(fit_hier2)) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  select(model, elpd_diff, se_diff) |>
  tt()

#' We get warnings about Pareto k's > 0.7 in PSIS-LOO, but as the
#' models are similar, and the difference is small, we can be
#' relatively confident that the more complex model is not better. In
#' this case, the likely reason is that the data do not have enough
#' information to learn about the interactions and adding them just
#' increases the posterior uncertainty.
#' 

#' <br />
#' 
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2017-2024, Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2017-2024, Aki Vehtari, licensed under CC-BY-NC 4.0.
