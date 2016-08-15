# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

library(ggplot2)
library(gridExtra)
library(tidyr)

# Hierarchical model for SAT-example data (BDA3, p. 102)

y <- c(28,8,-3,7,-1,1,18,12)
s <- c(15,10,16,11,9,11,10,18)

# Plot data, use normpdf to plot both the y_j and sigma_j
x <- seq(-40, 60, length.out = 500)

df_sep <- mapply(function(y, s, x) dnorm(x, y, s), y, s, MoreArgs = list(x = x)) %>%
  as.data.frame() %>% setNames(LETTERS[1:8]) %>% cbind(x) %>% gather(school, p, -x)

labs1 <- c('Other Schools', 'School A')
plot_sep <- ggplot(data = df_sep) +
  geom_line(aes(x = x, y = p, color = (school=='A'), group = school)) +
  labs(x = 'Treatment effect', y = '', title = 'Separate model', color = '') +
  scale_y_continuous(breaks = NULL) +
  scale_color_manual(values = c('blue','red'), labels = labs1) +
  theme(legend.background = element_blank(), legend.position = c(0.8,0.9))
plot_sep

df_pool <- data.frame(x = x, p = dnorm(x, sum(y/s^2)/sum(1/s^2), sqrt(1/sum(1/s^2))))

plot_pool <- ggplot(data = df_pool) +
  geom_line(aes(x = x, y = p, color = '1')) +
  labs(x = 'Treatment effect', y = '', title = 'Pooled model', color = '') +
  scale_y_continuous(breaks = NULL) +
  scale_color_manual(values = 'red', labels = 'All schools') +
  theme(legend.background = element_blank(), legend.position = c(0.7,0.9))

# Load the pre-computed results for the hierarchical model
# Replace this with your own code in the related exercise
load('demo5_2.RData')

df_hier <- as.data.frame(t(pxm)) %>% setNames(LETTERS[1:8]) %>%
  cbind(x) %>% gather(school, p, -x)

plot_hier <- ggplot(data = df_hier) +
  geom_line(aes(x = x, y = p, color = (school=='A'), group = school)) +
  labs(x = 'Treatment effect', y = '', title = 'Hierarchical model', color = '') +
  scale_y_continuous(breaks = NULL) +
  scale_color_manual(values = c('blue','red'), labels = labs1) +
  theme(legend.background = element_blank(), legend.position = c(0.8,0.9))

# Plot separate, pooled, and hierarchical model
grid.arrange(plot_sep, plot_pool, plot_hier)

# Various marginal and conditional posterior summaries

df_margpost = data.frame(x = t(tt), p = t(tp))

title1 <- 'Marginal posterior density p(tau|y)'
plot_margpost <-
  ggplot(data = df_margpost) +
  geom_line(aes(x = x, y = p)) +
  labs(x = expression(tau), y = 'p(tau|y)', title = title1) +
  scale_y_continuous(breaks = NULL)

df_condmeans <- as.data.frame(t(tm)) %>% setNames(LETTERS[1:8]) %>%
  cbind(x = t(tt)) %>% gather(school, p, -x)

yl <- c(-50, 50)
title2 <- 'Conditional posterior means of effects E[theta_j|tau,y]'
plot_condmeans <- ggplot(data = df_condmeans) +
  geom_line(aes(x = x, y = p, color = (school=='A'), group = school)) +
  coord_cartesian(ylim = yl) +
  labs(x = expression(tau), y = 'E[theta_j|tau,y)', title = title2, color = '') +
  scale_color_manual(values = c('blue','red'), labels = labs1) +
  theme(legend.background = element_blank(), legend.position = c(0.8,0.9))

df_condsds <- as.data.frame(t(tsd)) %>% setNames(LETTERS[1:8]) %>%
  cbind(x = t(tt)) %>% gather(school, p, -x)

title3 <- 'Conditional posterior standard deviations of effects sd[theta_j|tau,y]'
plot_condsds <- ggplot(data = df_condsds) +
  geom_line(aes(x = x, y = p, color = (school=='A'), group = school)) +
  coord_cartesian(ylim = yl) +
  labs(x = expression(tau), y = 'sd[theta_j|tau,y)', title = title3, color = '') +
  scale_color_manual(values = c('blue','red'), labels = labs1) +
  theme(legend.background = element_blank(), legend.position = c(0.8,0.9))

# Plot the posterior summaries
grid.arrange(plot_margpost, plot_condmeans, plot_condsds)

