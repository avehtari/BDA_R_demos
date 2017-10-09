# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# Metropolis algorithm + PSRF demonstration

library(ggplot2)
library(gridExtra)
library(tidyr)
library(devtools)
# install_github("dgrtwo/gganimate",
#                ref = '26ec501b78e0853134a3fe50a04364aef13d5f6c')
library(gganimate)
library(MASS)
library(plyr)
library(coda)

# Parameters of a Normal distribution used as a toy target distribution
y1 <- 0
y2 <- 0
r <- 0.8
S <- diag(2)
S[1, 2] <- r
S[2, 1] <- r

# draw samples from the toy distribution to visualize 90% HPD
# interval with ggplot's stat_ellipse()
dft <- data.frame(mvrnorm(100000, c(0, 0), S))

# load pre-run Metropolis chains
# note that proposal distribution was intentionally selected to be
# slightly too small, to better illustrate convergence diagonstics
# Since, implementation of the Metropolis algorithm is on of the
# exercises, we load here pre-computed chains and PSRF-values

# `tts' contains samples, `p1' and 'p2' contain PSRF values for t1
# and t2 using 50% warm-up, `pp1' and 'pp2' contain PSRF values for
# t1 and t2 using 10% warm-up, PSRF-values have been computed for
# each time-step.

data_path <- 'demos_ch11/demo11_4.RData'
load(data_path)

# transform the first s1 rows of the
# data into a 'tidy' format for plotting
s1 <- 50
dfs1 <- data.frame(tts[1:s1, 1, ]) %>% gather(chain, th1) %>%
  within({th2 <- c(tts[1:s1, 2, ])         # xl and yl specify the previous
          th2l <- c(th2[1], th2[-length(th2)])   # observation in the chain for
          th1l <- c(th1[1], th1[-length(th1)])}) # plotting

# Fix the incorrect lagged values, the lagged value of the
# first sample in the chain (for plotting) is  the value
# itself (instead of the last value of the previous chain)
sind <- 0:9*s1+1
dfs1[sind, c('th1l','th2l')] <- dfs1[sind, c('th1','th2')]

# another data frame with observations 501-10000
inds2 <- 501:10000
dfs2 <- data.frame(iteration = inds2, tts[inds2, 1, ]) %>%
  gather(chain, theta1, -iteration) %>%
  within(theta2 <- c(tts[inds2, 2, ])) %>%
  gather(var, val, -iteration, -chain)

# third data frame with PSRF values
indsp <- seq(10, length(p1), 10)
dfp <- data.frame(iteration = indsp,
                  htheta1 = p1[indsp], ttheta1 = pp1[indsp],
                  htheta2 = p2[indsp], ttheta2 = pp2[indsp]) %>%
  gather(grp,psrf,-iteration) %>% separate(grp, c('warm','var'), 1)
# h refers to 50% warm-up and t to 10% warm-up

# construct a 2d-plot of the 50 first iterations of the chains
chains1 <- ggplot(data = dfs1) +
  geom_segment(aes(x = th1, xend = th1l, y = th2, yend = th2l, color = chain,
                   frame = rep(1:s1, 10), group = chain, cumulative = T)) +
  geom_point(data = dfs1[sind, ], aes(x = th1, y = th2, color = chain)) +
  stat_ellipse(data = dft, aes(x = X1, y = X2), level = 0.9, color = 'black') +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) +
  labs(x = 'theta1', y = 'theta2') +
  scale_color_discrete(guide = FALSE)

# Animate s1 first iterations of the chains. Note, that
# at some points some of the chains seem to half for
# a moment. What really happens at that point is that
# they draw possibly a few points that are rejected
# (rejected points not shown) and thus the chain is not moving
gg_animate(chains1, interval = 0.2)

# Plot the result, no convergence yet
chains1 + labs(title = 'No convergence')

# Construct 1d-plots of the chains
chains21 <- ggplot(data = dfs1) +
  geom_line(aes(x = rep(1:s1, 10), y = th1, color = chain)) +
  labs(title = 'Convergence demonstration - not yet converged', x = '') +
  scale_color_discrete(guide = FALSE)

chains22 <- ggplot(data = dfs1) +
  geom_line(aes(x = rep(1:s1, 10), y = th2, color = chain)) +
  labs(x = 'iteration') +
  scale_color_discrete(guide = FALSE)

# Plot trends of the 50 first samples
grid.arrange(chains21, chains22)


# Plot trends from 500th sample to end
ggplot(data = dfs2) +
  geom_line(aes(iteration, val, color = chain)) +
  facet_grid(var~.) +
  labs(title = 'Convergence demonstration - visually converged', y = '') +
  scale_color_discrete(guide = FALSE)

# Plot PSRF with 50% warm-up and 10% warm-up
ggplot(data = dfp) +
  geom_line(aes(iteration, psrf, color = warm)) +
  facet_grid(var~.) +
  geom_hline(aes(yintercept = 1), linetype = 'dashed') +
  labs(title = 'Running PSRF with different warmp-up length', y = '') +
  scale_x_log10(breaks = 10^(0:4)) +
  scale_color_discrete(labels = c('PSRF(n/2:n)','PSRF(n/10:n)')) +
  theme(legend.position = 'bottom', legend.title = element_blank())

# Demonstrate how to compute PSRF/Rhat using coda
# Note that coda is using the older version of PSRF
# We need to form an object supported by coda package
ttsm <- mcmc.list(alply(tts,3,mcmc))
# PSRF per iteration
gelman.plot(ttsm)
# total PSRF
gelman.diag(ttsm)
