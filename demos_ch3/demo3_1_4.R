# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# Examples and illustrations for a normal model with unknown mean and
# variance (BDA3 section 3.2 on p. 64).
# Multivariate joint distribution, conditional distribution, marginal
# distribution, marginalization and posterior predictive distribution
# R versions of demos 3_1-3_4 combined into one demo
# like iPython notebook versions

library(ggplot2)
library(grid)
library(gridExtra)
library(tidyr)

# data
y <- c(93, 112, 122, 135, 122, 150, 118, 90, 124, 114)
# sufficient statistics
n <- length(y)
s2 <- var(y)
my <- mean(y)

# Factorize the joint posterior p(mu,sigma2|y) to p(sigma2|y)p(mu|sigma2,y)
# Sample from the joint posterior using this factorization

# helper function to sample from and evaluate
# scaled inverse chi-squared distribution
rsinvchisq <- function(n, nu, s2, ...) nu*s2 / rchisq(n , nu, ...)

dsinvchisq <- function(x, nu, s2){
  exp(log(nu/2)*nu/2 - lgamma(nu/2) + log(s2)/2*nu - log(x)*(nu/2+1) - (nu*s2/2)/x)
}

# sample 1000 random numbers from p(sigma2|y)
ns <- 1000
sigma2  <- rsinvchisq(ns, n-1, s2)

# sample from p(mu|sigma2,y)
mu <- my + sqrt(sigma2/n)*rnorm(length(sigma2))

# construct a variable sigma and
sigma <- sqrt(sigma2)

# sample from predictive distribution p(ynew|y) for each sample of (mu, sigma)
ynew <- rnorm(ns, mu, sigma)

# For mu, sigma and ynew compute the density in a grid
# ranges for the grids
t1l <- c(90, 150)
t2l <- c(10, 60)
nl <- c(50, 185)

t1 <- seq(t1l[1], t1l[2], length.out = ns)
t2 <- seq(t2l[1], t2l[2], length.out = ns)
xynew <- seq(nl[1], nl[2], length.out = ns)

# compute the exact marginal density of mu
# multiplication by 1./sqrt(s2/n) is due to the transformation of
# variable z=(x-mean(y))/sqrt(s2/n), see Gelman et al p. 24
pm <- dt((t1-my) / sqrt(s2/n), n-1) / sqrt(s2/n)

# estimate the marginal density using samples
# and ad hoc Gaussian kernel approximation
pmk <- density(mu, adjust = 2, n = ns, from = t1l[1], to = t1l[2])$y

# compute the exact marginal density of sigma
# the multiplication by 2*t2 is due to the transformation of
# variable z=t2^2, see Gelman et al p. 24
ps <- dsinvchisq(t2^2, n-1, s2) * 2*t2

# estimate the marginal density using samples
# and ad hoc Gaussian kernel approximation
psk <- density(sigma, n = ns, from = t2l[1], to = t2l[2])$y

# compute the exact predictive density
# multiplication by 1./sqrt(s2/n) is due to the transformation of variable
# see BDA3 p. 21
p_new <- dt((xynew-my) / sqrt(s2*(1+1/n)), n-1) / sqrt(s2*(1+1/n))

# combine grid points into another data frame
# with all pairwise combinations
dfj <- data.frame(t1 = rep(t1, each = length(t2)),
                  t2 = rep(t2, length(t1)))

# evaluate the joint density in a grid
# note that the following is not normalized, but for plotting
# contours it does not matter
dfj$z <- dsinvchisq(dfj$t2^2, n-1, s2) * 2*dfj$t2 * dnorm(dfj$t1, my, dfj$t2/sqrt(n))
# breaks for plotting the contours
cl <- seq(1e-5, max(dfj$z), length.out = 6)

## Visualise the joint density and marginal densities of the posterior
## of normal distribution with unknown mean and variance.

dfm <- data.frame(t1, Exact = pm, Empirical = pmk) %>% gather(grp, p, -t1)
# create a plot of the marginal density of mu
margmu <- ggplot(dfm) +
  geom_line(aes(t1, p, color = grp)) +
  coord_cartesian(xlim = t1l) +
  labs(title = 'Marginal of mu', x = '', y = '') +
  scale_y_continuous(breaks = NULL) +
  theme(legend.background = element_blank(),
        legend.position = c(0.75, 0.8),
        legend.title = element_blank())

dfs <- data.frame(t2, Exact = ps, Empirical = psk) %>% gather(grp, p, -t2)
# create a plot of the marginal density of sigma
margsig <- ggplot(dfs) +
  geom_line(aes(t2, p, color = grp)) +
  coord_cartesian(xlim = t2l) +
  coord_flip() +
  labs(title = 'Marginal of sigma', x = '', y = '') +
  scale_y_continuous(breaks = NULL) +
  theme(legend.background = element_blank(),
        legend.position = c(0.75, 0.8),
        legend.title = element_blank())

# create a plot of the joint density
joint1labs <- c('Samples','Exact contour')
joint1 <- ggplot() +
  geom_point(data = data.frame(mu,sigma), aes(mu, sigma, col = '1'), size = 0.1) +
  geom_contour(data = dfj, aes(t1, t2, z = z, col = '2'), breaks = cl) +
  coord_cartesian(xlim = t1l,ylim = t2l) +
  labs(title = 'Joint posterior', x = '', y = '') +
  scale_y_continuous(labels = NULL) +
  scale_x_continuous(labels = NULL) +
  scale_color_manual(values=c('blue', 'black'), labels = joint1labs) +
  guides(color = guide_legend(nrow  = 1, override.aes = list(
    shape = c(16, NA), linetype = c(0, 1), size = c(2, 1)))) +
  theme(legend.background = element_blank(),
        legend.position = c(0.5, 0.9),
        legend.title = element_blank())

# blank plot for combining the plots
bp <- grid.rect(gp = gpar(col = 'white'))

# combine the plots
grid.arrange(joint1, margsig, margmu, bp, nrow = 2)

## Visualise factored sampling and the corresponding
## marginal and conditional densities.

# data frame for the conditional of mu and marginal of sigma
dfc <- data.frame(mu = t1, marg = rep(sigma[1], length(t1)),
                  cond = sigma[1] + dnorm(t1 ,my, sqrt(sigma2[1]/n)) * 100) %>%
  gather(grp, p, marg, cond)

# legend labels for the following plot
joint2labs <- c('Exact contour plot', 'Sample from joint post.',
               'Cond. distribution of mu', 'Sample from the marg. of sigma')

# create another plot of the joint posterior
joint2 <- ggplot() +
  geom_contour(data = dfj, aes(t1, t2, z = z, col = '1'), breaks = cl) +
  geom_point(data = data.frame(m = mu[1], s = sigma[1]), aes(m , s, color = '2')) +
  geom_line(data = dfc, aes(mu, p, color = grp), linetype = 'dashed') +
  coord_cartesian(xlim = t1l,ylim = t2l) +
  labs(title = 'Joint posterior', x = '', y = '') +
  scale_x_continuous(labels = NULL) +
  scale_color_manual(values=c('blue', 'darkgreen','darkgreen','black'), labels = joint2labs) +
  guides(color = guide_legend(nrow  = 2, override.aes = list(
    shape = c(NA, 16, NA, NA), linetype = c(1, 0, 1, 1)))) +
  theme(legend.background = element_blank(),
        legend.position = c(0.5, 0.85),
        legend.title = element_blank())

# create another plot of the marginal density of sigma
margsig2 <- ggplot(data = data.frame(t2, ps)) +
  geom_line(aes(t2, ps), color = 'blue') +
  coord_cartesian(xlim = t2l) +
  coord_flip() +
  labs(title = 'Marginal of sigma', x = '', y = '') +
  scale_y_continuous(labels = NULL)

# combine the plots
grid.arrange(joint2, margsig2, ncol = 2)

## Visualise the marginal distribution of mu as a mixture of normals.

# calculate conditional pdfs for each sample
condpdfs <- sapply(t1, function(x) dnorm(x, my, sqrt(sigma2/n)))

# data frame of 25 first samples
dfm25 <- data.frame(t1, t(condpdfs[1:25,])) %>% gather(grp, p, -t1)

# create a plot of some of them
condmu <- ggplot(data = dfm25) +
  geom_line(aes(t1, p, group = grp), linetype = 'dashed') +
  labs(title = 'Conditional distribution of mu for first 25 samples', y = '', x = '') +
  scale_y_continuous(breaks = NULL)


# labels for the following plot
mulabs <- c('avg of sampled conds', 'exact marginal of mu')

dfsam <- data.frame(t1, colMeans(condpdfs), pm) %>% gather(grp,p,-t1)
# create a plot of their mean
meanmu <- ggplot(data = dfsam) +
  geom_line(aes(t1, p, size = grp, color = grp)) +
  labs(y = '', x = 'mu') +
  scale_y_continuous(breaks = NULL) +
  scale_size_manual(values = c(2, 0.8), labels = mulabs) +
  scale_color_manual(values = c('orange', 'black'), labels = mulabs) +
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_blank(),
        legend.title = element_blank())

# combine the plots
grid.arrange(condmu, meanmu, ncol = 1)

## Visualise sampling from the posterior predictive distribution.

# calculate each predictive pdf with given mu and sigma
ynewdists <- sapply(xynew, function(x) dnorm(x, mu, sigma))

# data frame of dirst draw from sample the predictive along
# with the exact value for plotting
dfp <- data.frame(xynew, draw = ynewdists[1,], exact = p_new)

# data frame for plotting the samples
dfy <- data.frame(ynew, p = 0.02*max(ynewdists))

# legend labels for the plots
pred1labs <- c('Sample from the predictive distribution', 'Predictive distribution given the posterior sample')
pred2labs <- c('Samples from the predictive distribution', 'Exact predictive distribution')

# create yet another plot of the joint posterior with a draw
# from the posterior predictive distribution, highlight the first sample
# create a plot of the joint density
joint3labs <- c('Samples', 'Exact contour')
joint3 <- ggplot() +
  geom_point(data = data.frame(mu, sigma), aes(mu, sigma, col = '1'), size = 0.1) +
  geom_contour(data = dfj, aes(t1, t2, z = z, col = '2'), breaks = cl) +
  geom_point(data = data.frame(x = mu[1], y = sigma[1]), aes(x, y), color = 'red') +
  coord_cartesian(xlim = t1l,ylim = t2l) +
  labs(title = 'Joint posterior', x = 'mu', y = 'sigma') +
  scale_color_manual(values=c('blue', 'black'), labels = joint3labs) +
  guides(color = guide_legend(nrow  = 1, override.aes = list(
    shape = c(16, NA), linetype = c(0, 1), size = c(2, 1)))) +
  theme(legend.background = element_blank(),
        legend.position=c(0.5 ,0.9),
        legend.title = element_blank())

# create a plot of the predicitive distribution and the respective sample
pred1 <- ggplot() +
  geom_line(data = dfp, aes(xynew, draw, color = '2')) +
  geom_point(data = dfy, aes(ynew[1], p, color = '1')) +
  coord_cartesian(xlim = nl, ylim = c(0,0.04)) +
  labs(title = 'Posterior predictive distribution', x = 'ytilde', y = '') +
  scale_y_continuous(breaks = NULL) +
  scale_color_manual(values = c('red', 'blue'), labels = pred1labs) +
  guides(color = guide_legend(nrow = 2, override.aes = list(
    linetype = c(0, 1), shape=c(16, NA), labels = pred1labs))) +
  theme(legend.background = element_blank(),
        legend.position = c(0.5 ,0.9),
        legend.title = element_blank())

# create a plot for all ynews
pred2 <- ggplot() +
  geom_line(data = dfp, aes(xynew, draw, color = '2')) +
  geom_point(data = dfy, aes(ynew, p, color = '1'), alpha = 0.05) +
  coord_cartesian(xlim = nl, ylim = c(0,0.04)) +
  labs(x = 'ytilde', y = '') +
  scale_y_continuous(breaks=NULL) +
  scale_color_manual(values=c('darkgreen','blue'),labels=pred2labs) +
  guides(color = guide_legend(nrow = 2, override.aes=list(
    linetype = c(0, 1), shape = c(16, NA), alpha = c(1, 1) ,labels = pred2labs))) +
  theme(legend.background = element_blank(),
        legend.position = c(0.5 ,0.9),
        legend.title = element_blank())

# combine the plots
grid.arrange(joint3, pred1, bp, pred2, nrow = 2)

