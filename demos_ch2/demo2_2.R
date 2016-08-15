# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# Illustrate the effect of prior.
# Comparison of posterior distributions with different
# parameter values for the beta prior distribution.

# ggplot2 is used for plotting, tidyr for manipulating data frames
library(ggplot2)
library(tidyr)

# Prior, posterior and data:
# observed data, 437 girls and 543 boys
a <- 437
b <- 543

# Evaluate densities at evenly spaced points between 0.375 and 0.525
df1 <- data.frame(x = seq(0.375, 0.525, 0.001))

# Posterior with Beta(1,1), ie. uniform prior
df1$pu <- dbeta(df1$x, a+1, b+1)

# 3 different choices for priors, Beta(0.485*2,(1-0.485)*2),
# Beta(0.485*20,(1-0.485)*20) and Beta(0.485*200,(1-0.485)*200)
n <- c(2, 20, 200) # prior counts
apr <- 0.485 # prior ratio of success

# helperf returns for given number of prior observations, prior ratio
# of successes, number of observed successes and failures and a data
# frame with values of x, a new data frame with prior and posterior
# values evaluated at points x.
helperf <- function(n, apr, a, b, df)
  cbind(df, pr = dbeta(df$x, n*apr, n*(1-apr)), po = dbeta(df$x, n*apr + a, n*(1-apr) + b), n = n)

# lapply function over prior counts n and gather results into key-value pairs.
df2 <- lapply(n, helperf, apr, a, b, df1) %>% do.call(rbind, args = .) %>%
  gather(grp, p, -c(x, n), factor_key = T)

# add correct labels for plotting
df2$title <- factor(paste0('alpha/(alpha+beta)=0.485, alpha+beta=',df2$n))
levels(df2$grp) <- c('Posterior with unif prior', 'Prior', 'Posterior')

# plot distributions
ggplot(data = df2) +
  geom_line(aes(x, p, color = grp)) +
  geom_vline(xintercept = 0.485, linetype = 'dotted') +
  facet_wrap(~title, ncol = 1) +
  labs(x = '', y = '') +
  scale_y_continuous(breaks = NULL) +
  theme(legend.position = 'bottom', legend.title = element_blank())

