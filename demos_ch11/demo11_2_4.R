#' ---
#' title: "Bayesian data analysis demos 11.2, 11.3, and 11.4"
#' author: "Aki Vehtari, Markus Paasiniemi"
#' date: "`r format(Sys.Date())`"
#' ---

#' ## Metropolis algorithm + Rhat (PSRF) demonstration

#' ggplot2 and gridExtra are used for plotting, tidyr for manipulating
#' data frames
#+ setup, message=FALSE, error=FALSE, warning=FALSE
library(ggplot2)
theme_set(theme_minimal())
library(gridExtra)
library(tidyr)
#library(devtools)
#install_github("dgrtwo/gganimate")
library(gganimate)
library(MASS)
library(plyr)
library(rstan)
library(rprojroot)
root<-has_dirname("BDA_R_demos")$make_fix_file()

#' Parameters of a Normal distribution used as a toy target distribution
y1 <- 0
y2 <- 0
r <- 0.8
S <- diag(2)
S[1, 2] <- r
S[2, 1] <- r

#' Draw samples from the toy distribution to visualize 90% HPD
#' interval with ggplot's stat_ellipse()
dft <- data.frame(mvrnorm(100000, c(0, 0), S))

#' ### Load pre-run Metropolis chains.
#' Since, implementation of the Metropolis algorithm is one of the
#' exercises, we load here pre-computed chains and Rhat-values
#' 
#' The proposal distribution was intentionally selected to be slightly
#' too small, to better illustrate convergence diagonstics
#' 
#' `tts' contains samples, `p1' and 'p2' contain Rhat values for t1
#' and t2 using 50% warm-up, `pp1' and 'pp2' contain Rhat values for
#' t1 and t2 using 10% warm-up, Rhat-values have been computed for
#' each time-step.
load(root("demos_ch11","demo11_4.RData"))

#' Transform the first s1 rows of the
#' data into a 'tidy' format for plotting
s1 <- 50
dfs1 <- data.frame(iter = 1:s1, tts[1:s1, 1, ]) %>%
  gather(chain, th1, -iter) %>%
  within({th2 <- c(tts[1:s1, 2, ])         # xl and yl specify the previous
          th2l <- c(th2[1], th2[-length(th2)])   # observation in the chain for
          th1l <- c(th1[1], th1[-length(th1)])}) # plotting

#' Fix the incorrect lagged values, the lagged value of the
#' first sample in the chain (for plotting) is  the value
#' itself (instead of the last value of the previous chain)
sind <- 0:9*s1+1
dfs1[sind, c('th1l','th2l')] <- dfs1[sind, c('th1','th2')]

#' Another data frame with observations 501-10000
inds2 <- 501:10000
dfs2 <- data.frame(iter = inds2, tts[inds2, 1, ]) %>%
  gather(chain, theta1, -iter) %>%
  within(theta2 <- c(tts[inds2, 2, ])) %>%
  gather(var, val, -iter, -chain)

#' Third data frame with Rhat values
indsp <- seq(10, length(p1), 10)
dfp <- data.frame(iter = indsp,
                  theta1 = p1[indsp],
                  theta2 = p2[indsp]) %>%
  gather(var,rhat,-iter) 

#' Construct a 2d-plot of the 50 first iterations of the chains
frame = rep(1:s1, 10)
chains1 <- ggplot(data = dfs1) +
  geom_segment(aes(x = th1, xend = th1l, y = th2, yend = th2l, color = chain,
                   group = chain)) +
  geom_point(data = dfs1[sind, ], aes(x = th1, y = th2, color = chain)) +
  stat_ellipse(data = dft, aes(x = X1, y = X2), level = 0.9, color = 'black') +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) +
  labs(x = 'theta1', y = 'theta2') +
  scale_color_discrete(guide = FALSE)

#' Animate s1 first iterations of the chains. 
#' At some points some of the chains seem to halt for
#ä a moment. What really happens at that point is that
#' they draw possibly a few points that are rejected
#' (rejected points not shown) and thus the chain is not moving
#+ Metropolis
animate(chains1 + 
          transition_reveal(id=chain, along=iter) + 
          shadow_trail(1/s1))

#' Plot the result, no convergence yet
chains1 + labs(title = 'No convergence')

#' Construct 1d-plots of the chains
chains21 <- ggplot(data = dfs1) +
  geom_line(aes(x = rep(1:s1, 10), y = th1, color = chain)) +
  labs(title = 'Convergence demonstration - not yet converged', x = '') +
  scale_color_discrete(guide = FALSE)
chains22 <- ggplot(data = dfs1) +
  geom_line(aes(x = rep(1:s1, 10), y = th2, color = chain)) +
  labs(x = 'iteration') +
  scale_color_discrete(guide = FALSE)

#' Plot trends of the 50 first samples
grid.arrange(chains21, chains22)

#' Plot trends from 500th sample to end
ggplot(data = dfs2) +
  geom_line(aes(iter, val, color = chain)) +
  facet_grid(var~.) +
  labs(title = 'Convergence demonstration - visually converged', y = '') +
  scale_color_discrete(guide = FALSE)

#' Plot Rhat with 50% warm-up and increasing number of iterations
ggplot(data = dfp) +
  geom_line(aes(iter, rhat, color = var)) +
  geom_hline(aes(yintercept = 1), linetype = 'dashed') +
  labs(title = 'Running Rhat with different warmp-up length', y = '') +
  scale_x_log10(breaks = 10^(0:4)) +
  scale_color_discrete(labels = c('theta1','theta2')) +
  theme(legend.position = 'bottom', legend.title = element_blank())

#' Demonstrate how to compute Rhat (PSRF) using RStan monitor.
#' We need to reorder our array
samp <- aperm(tts, c(1, 3, 2))
monitor(samp, probs = c(0.25, 0.5, 0.75), digits_summary = 2)
