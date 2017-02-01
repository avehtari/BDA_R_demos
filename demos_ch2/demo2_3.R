# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# Probability of a girl birth given placenta previa (BDA3 p. 37).
# Simulate samples from Beta(438,544), draw a histogram with
# quantiles, and do the same for a transformed variable.

library(ggplot2)
library(tidyr)

# Sample from posterior Beta(438,544)
# Get all samples at once and store them in vector 'th'
a <- 438
b <- 544

th <- rbeta(10000, a, b)
# phi = (1-theta)/theta
phi <- (1 - th) / th

# merge the data into one data frame for plotting
df1 <- data.frame(phi,th) %>% gather()

# compute 2.5% and 97.5% quantile approximation using samples
quantiles <- c(0.025, 0.975)
thq <- quantile(th, quantiles)
phiq <- quantile(phi, quantiles)

# merge quantiles into one data frame for plotting
df2 <- data.frame(phi = phiq, th = thq) %>% gather()

# Histogram plots with 30 bins for theta and phi
ggplot() +
  geom_histogram(data = df1, aes(value), bins = 30) +
  geom_vline(data = df2, aes(xintercept = value), linetype = 'dotted') +
  facet_wrap(~key, ncol = 1, scales = 'free_x')  +
  labs(x = '', y = '') +
  scale_y_continuous(breaks = NULL)

