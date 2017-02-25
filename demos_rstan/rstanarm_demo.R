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

# The following models do not equal the models
# at rstan_demo.R exactly, but rather serve as
# examples of how to implement similar models
# with rstanarm


# Bernoulli model
d_bern <- list(y = c(0, 1, 0, 0, 1, 1, 1, 0, 1, 0))

# With a uniform prior (beta(1,1)).
# This is achieved by setting the prior to NULL,
# which is not recommended in general.
# y ~ 1 means y depends only on the intercept term
fit_bern <- stan_glm(y ~ 1, family = binomial(),
                     data = d_bern, prior_intercept = NULL)

# One convenient way to examine and diagnose the
# fitted model is to call shinystan as follows:
launch_shinystan(fit_bern)

# to see the parameter values on the ouput space, do the
# inverse logistic transformation (plogis in R) on them
# intercept
coef(fit_bern)
# probability of success
plogis(coef(fit_bern))

#intercept
stan_hist(fit_bern)
# probability of success
theta <- plogis(extract(fit_bern$stanfit)$alpha)
ggplot() + geom_histogram(aes(theta), bins = 50, fill = 'darkblue', color = 'black')

# Binomial model with a roughly uniform prior for
# the probability of success.
# The prior is specified in the 'latent space'. The
# actual probability of success, theta = plogis(alpha),
# where plogis is the inverse of the logistic function.
#
# Visualize the prior by drawing samples from it
prior_mean <- 0
prior_sd <- 1.5
prior_intercept <- normal(location = prior_mean, scale = prior_sd)
prior_samples <- plogis(rnorm(20000, prior_mean, prior_sd))
ggplot() + geom_histogram(aes(prior_samples), bins = 25, fill = 'darkblue', color = 'black')

d_bin <- list(N = c(5,5), y = c(4,3))
fit_bin <- stan_glm(y/N ~ 1, family = binomial(), data = d_bin,
                     prior_intercept = prior_intercept, weights = N)
#launch_shinystan(fit_bern)
plogis(coef(fit_bin))
ggplot() + geom_histogram(aes(x = plogis(extract(fit_bin$stanfit)$alpha)),
                          bins = 50, fill = 'darkblue', color = 'black') +
  labs(x = 'probability of success', y = '') + scale_y_continuous(breaks = NULL)

# Re-run the model with a new data dataset.
d_bin <- list(N = c(5,5), y = c(4,5))
fit_bin <- update(fit_bin, data = d_bin)
#launch_shinystan(fit_bern)
plogis(coef(fit_bin))
ggplot() + geom_histogram(aes(x = plogis(extract(fit_bin$stanfit)$alpha)),
                          bins = 50, fill = 'darkblue', color = 'black') +
  labs(x = 'probability of success', y = '') + scale_y_continuous(breaks = NULL)

# comparison of two groups with Binomial
# grp2 is a dummy variable that captures the
# differece of the intercepts in the first and
# the second group
d_bin2 <- list(N = c(674, 680), y = c(39,22), grp2 = c(0,1))
fit_bin2 <- stan_glm(y/N ~ grp2, family = binomial(), data = d_bin2,
                     prior_intercept = NULL, prior = NULL, weights = N)
samples_bin2 <- extract(fit_bin2$stanfit)
theta1 <- plogis(samples_bin2$alpha)
theta2 <- plogis(samples_bin2$alpha + samples_bin2$beta)
oddsratio <- (theta2/(1-theta2))/(theta1/(1-theta1))
ggplot() + geom_histogram(aes(oddsratio), bins = 50, fill = 'darkblue', color = 'black') +
  labs(y = '') + scale_y_continuous(breaks = NULL)


# Linear model
d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_lin <-list(year = d_kilpis$year,
             temp = d_kilpis[,5])

# y ~ x means y depends on the intercept and x
fit_lin <- stan_glm(temp ~ year, data = d_lin, family = gaussian())

#launch_shinystan(fit_lin)
samples_lin <- rstan::extract(fit_lin$stanfit, permuted = T)
mean(samples_lin$beta>0) # probability that beta > 0
mu_samples <- tcrossprod(cbind(1, d_lin$year), cbind(samples_lin$alpha,samples_lin$beta))

mu <- apply(mu_samples, 1, quantile, c(0.05, 0.5, 0.95)) %>%
  t() %>% data.frame(x = d_lin$year, .) %>% gather(pct, y, -x)
pfit <- ggplot() +
  geom_point(aes(year, temp), data = data.frame(d_lin), size = 0.5) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(x = '', y = 'Summer temp. @Kilpisjärvi') +
  guides(linetype = F) +
  theme_bw()
phist <- stan_hist(fit_lin, pars = c('beta','sigma'), bins = 50) + ggtitle('parameters')
grid.arrange(pfit, phist)

# psis-loo
loo1 <- loo(fit_lin)

# prediction for a new data point
predict(fit_lin, newdata = data.frame(year = 2016), se.fit = T)
# or sample from the posterior predictive distribution and
# plot the histogram
post_samples <- posterior_predict(fit_lin, newdata = data.frame(year = 2016))
ggplot(data = data.frame(ypred = post_samples)) +
  geom_histogram(aes(ypred), bins = 50, fill = 'darkblue', color = 'black') +
  labs(y = '', x = 'avg-temperature prediction for the summer 2016') +
  scale_y_continuous(breaks = NULL)

# Currently, rstanarm does not support student-t likelihood

# comparison of k groups

d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_grp <- list(month = rep(6:8, nrow(d_kilpis)),
              temp = c(t(d_kilpis[,2:4])))

# weakly informative prior for the common mean
prior_intercept <- normal(10, 10)
# again, to use no (= uniform) prior, prior_intercept
# could be set to NULL

# y ~ 1 + (1 | x) means y depends on common intercept and
# group speficific intercepts (grouping determined by x)
fit_grp <- stan_lmer(temp ~ 1 + (1 | month), data = d_grp,
                     prior_intercept = prior_intercept)
#launch_shinystan(fit_grp)
# average temperature and monthly
# deviations from the mean
stan_hist(fit_grp, bins = 50)

# A boxplot like the one in rstan_demo.R
# can be obtained as follows:
temps <- (as.matrix(fit_grp)[,1] + as.matrix(fit_grp)[, 2:4]) %>%
  as.data.frame() %>% setNames(6:8) %>% gather(month, temp)
qplot(month, temp, data = temps, geom = 'boxplot')
# or a  similar plot:
# stan_plot(fit_grp)

# probabilities that june is hotter than july, june is hotter than august
# and july is hotter than august:
combn(unique(temps$month), 2, function(months, data) {
  mean(subset(data, month == months[1])$temp > subset(data, month == months[2])$temp)
}, data = temps) %>% setNames(c('6>7', '6>8', '7>8'))

