# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# rstanarm demos

library(tidyr)
library(rstanarm)
library(rstan)
library(loo)
library(shinystan)
library(ggplot2)
library(gridExtra)

# the following models do not equal the models
# at rstan_dem.R exactly, but rather serve as
# examples of how to implement similar models
# with rstanarm
# # # #
# Bernoulli model with (improper) uniform prior
d_bern <- list(y = c(0, 1, 0, 0, 1, 1, 1, 0, 1, 0))
family <- binomial()
# y ~ 1 means y depends only on the intercept term
fit_bern <- stan_glm(y ~ 1, family = family, data = d_bern,
                     prior_intercept = NULL)

# one convenient way to examine and diagnose the
# fitted model is to call shinystan as follows:
# (to see parameter values on ouput space, select
#  inv_logit-transformation in the diagnose-tab)
launch_shinystan(fit_bern)

# Otherwise, posterior samples can be extracted
# from the stanfit object
samples_bern <- extract(fit_bern$stanfit)

# when using other glms than basic linear model
# (like bernoulli in this case), transform the
# coefficients from the latent space to the output
# space with family$linkinv()
theta <- family$linkinv(samples_bern$alpha)
qplot(x = theta, geom = 'histogram', bins = 50)

# # # #
# Bernoulli model with weakly informative prior
d_bern <- list(y = c(0, 1, 0, 0, 1, 1, 1, 0, 1, 0))
family <- binomial()
prior_intercept = normal(location = 0, scale = 2.5)
# scale 2.5 means, that approximately 95%
# of prior mass for theta is between values
c(family$linkinv(0-2*2.5), family$linkinv(0+2*2.5))
# (the prior mean and mode being at 0.5).

fit_bern <- stan_glm(y ~ 1, family = family, data = d_bern,
                     prior_intercept = prior_intercept)
launch_shinystan(fit_bern)

# # # #
# Same model with binomial (as opposed to bernoulli) data
family <- binomial()
d_bin <- list(N = c(10,6), y = c(8,3))
prior_intercept = normal(location = 0, scale = 2.5)
fit_bern <- stan_glm(y/N ~ 1, family = family, data = d_bin,
                     prior_intercept = prior_intercept, weights = N)
launch_shinystan(fit_bern)

# # # #
# Gaussian linear model with (improper) uniform priors
d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_lin <- list(x = rep(d_kilpis$year, each = 3), xpred = 2016,
              y = c(t(d_kilpis[,2:4])), N = 3*nrow(d_kilpis))

# y ~ x means y depends on the intercept and x
fit_lin <- stan_glm(y ~ x, data = d_lin, prior = NULL,
                    prior_intercept = NULL)
launch_shinystan(fit_lin)
# probability that beta > 0
mean(extract(fit_lin$stanfit)$beta>0)


# psis-loo
loo1 <- loo(fit_lin)

# # # #
# Gaussian linear model with adjustable priors
d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_lin <- list(x = rep(d_kilpis$year, each = 3), xpred = 2016,
              y = c(t(d_kilpis[,2:4])), N = 3*nrow(d_kilpis))

pmualpha <- mean(d_lin$y) # centered
psalpha <- (14-4)/6       # avg temp between 4-14
pmubeta <- 0              # a priori incr. and decr. as likely
psbeta <- (.1--.1)/6      # avg temp prob does does not
prior_intercept <- normal(pmualpha, psalpha)
prior <- normal(pmubeta, psbeta)
fit_lin <- stan_glm(y ~ x, data = d_lin,  prior = prior,
                    prior_intercept = prior_intercept)
launch_shinystan(fit_lin)
# probability that beta > 0
mean(extract(fit_lin$stanfit)$beta>0)

# # # #
# comparison of k groups with common variance (ANOVA)
d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_grp <-list(x = rep(1:3, nrow(d_kilpis)),
             y = c(t(d_kilpis[,2:4])))

# weakly informative prior for the common mean
prior_intercept <- normal(10, 10)

# y ~ (1 | x) means y depends on common intercept and
# group speficific intercepts (grouping determined by x)
fit_grp <- stan_lmer(y ~ (1 | x), data = d_grp,
                     prior_intercept = prior_intercept)

launch_shinystan(fit_lin)

# plot(fit_grp) plots the common mean and the
# deviations from the common mean. A boxplot
# like the one in rstan_demo.R can be obtained
# as. follows:
mu <- (as.matrix(fit_grp)[,1] + as.matrix(fit_grp)[, 2:4]) %>%
   as.data.frame() %>% setNames(6:8) %>% gather(month, temp)
qplot(month, temp, data = mu, geom = 'boxplot')

# probabilities that june is hotter than july, june is hotter than august
# and july is hotter than august:
combn(unique(mu$month), 2, function(a, df)
  mean(subset(df, month == a[1])$temp > subset(df, month == a[2])$temp), df = mu)

