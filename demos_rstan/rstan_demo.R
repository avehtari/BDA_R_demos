# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# Stan demos

library(tidyr) #
library(rstan) # version >= 2.11
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(ggplot2)
library(gridExtra)

# Note that the stan-models are stored in separate .stan-files.


# Bernoulli model
d_bern <- list(N = 10, y = c(0, 1, 0, 0, 1, 1, 1, 0, 1, 0))

# With a beta(1,1) (uniform) prior
(fit_bern <- stan(file = 'bern.stan', data = d_bern))
stan_hist(fit_bern, bins = 50)
# or extract the samples for plotting manually:
# samples_bern <- extract(fit_bern, permuted = T)
# hist(samples_bern$theta)

# Binomial model with a roughly uniform prior for
# the probability of success.
# The prior is specified in the 'latent space'. The
# actual probability of success, theta = plogis(alpha),
# where plogis is the inverse of the logistic function.
#
# Visualize the prior by drawing samples from it
prior_samples <- plogis(rnorm(20000, 0, 1.5))
ggplot() + geom_histogram(aes(prior_samples), bins = 50, fill = 'darkblue', color = 'black')

d_bin <- list(N = 10, y = 7)
(fit_bin <- stan(file = 'binom.stan', data = d_bin))
stan_hist(fit_bin, pars = 'theta', bins = 50)

# Re-run the model with a new data dataset.
d_bin <- list(N = 10, y = 9)
(fit_bin <- stan(file = 'binom.stan', data = d_bin))
stan_hist(fit_bin, pars = 'theta', bins = 50)

# Comparison of two groups with Binomial
d_bin2 <- list(N1 = 674, y1 = 39, N2 = 680, y2 = 22)
(fit_bin2 <- stan(file = 'binom2.stan', data = d_bin2))
stan_hist(fit_bin2, pars = 'oddsratio', bins = 50)+geom_vline(xintercept = 1)

# Linear model example with Kilpisjärvi data
d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_lin <-list(N = nrow(d_kilpis),
             x = d_kilpis$year,
             xpred = 2016,
             y = d_kilpis[,5])

# Plot the data
ggplot() +
  geom_point(aes(x, y), data = data.frame(d_lin), size = 0.5) +
  labs(y = 'Summer temp. @Kilpisjärvi', x= "Year") +
  guides(linetype = F) +
  theme_bw()

# create another list with data and priors
d_lin_priors <- c(list(
    pmualpha = mean(unlist(d_kilpis[,5])), # centered
    psalpha = (14-4)/6, # avg temp between 4-14
    pmubeta = 0, # a priori incr. and decr. as likely
    psbeta = (.1--.1)/6), # avg temp prob does does not incr. more than a degree per 10 years
  d_lin)

# Gaussian linear model
# with adjustable priors
fit_lin <- stan(file = 'lin.stan', data = d_lin_priors)

# with standardized data
# this is alternative to above
#fit_lin <- stan(file = 'lin_std.stan', data = d_lin)

# Linear student-t model
fit_lin_t <- stan(file = 'lin_t.stan', data = d_lin)

samples_lin_t <- rstan::extract(fit_lin_t, permuted = T)
mean(samples_lin_t$beta>0) # probability that beta > 0
mu <- apply(samples_lin_t$mu, 2, quantile, c(0.05, 0.5, 0.95)) %>%
  t() %>% data.frame(x = d_lin$x, .)  %>% gather(pct, y, -x)

pfit <- ggplot() +
  geom_point(aes(x, y), data = data.frame(d_lin), size = 0.5) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(y = 'Summer temp. @Kilpisjärvi', x= "Year") +
  guides(linetype = F) +
  theme_bw()
pars <- intersect(names(samples_lin_t), c('beta','sigma','nu','ypred'))
phist <- stan_hist(fit_lin_t, pars = pars, bins = 50)
grid.arrange(pfit, phist, nrow = 2)

# psis-loo
# For the following three lines to execute, the log-likelihood
# needs to be evaluated in the stan code. For an example, see lin.stan
if('log_lik' %in% samples_lin_t) {
  log_lik <- extract_log_lik(fit_lin, parameter_name = 'log_lik')
  loo_lin <- loo(log_lik)
  log_lik_t <- extract_log_lik(fit_lin_t, parameter_name = 'log_lik')
  loo_lin_t <- loo(log_lik_t)
  compare(loo_lin,loo_lin_t)
}


# Comparison of k groups
d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_grp <-list(N = 3*nrow(d_kilpis),
             K = 3,
             x = rep(1:3, nrow(d_kilpis)),
             y = c(t(d_kilpis[,2:4])))

# common variance (ANOVA)
fit_grp <- stan(file = 'grp_aov.stan', data = d_grp)
fit_grp2 <- stan(file = 'grp_aov2.stan', data = d_grp)

# common variance and hierarchical prior for mean
# results do not differ much from the previous, because there is only
# few groups and quite much data per group, but this works as an example anyway
fit_grp <- stan(file = 'grp_prior_mean.stan', data = d_grp)

# unequal variance and hierarchical prior for mean and variance
fit_grp <- stan(file = 'grp_prior_mean_var.stan', data = d_grp)

# plot the results
samples_grp <- extract(fit_grp, permuted = T)

temps <- data.frame(samples_grp$mu) %>% setNames(6:8) %>% gather(month, temp)
qplot(month, temp, data = temps, geom = 'boxplot')

# probabilities that june is hotter than july, june is hotter than august
# and july is hotter than august:
combn(unique(temps$month), 2, function(months, data) {
  mean(subset(data, month == months[1])$temp > subset(data, month == months[2])$temp)
}, data = temps) %>% setNames(c('6>7', '6>8', '7>8'))
