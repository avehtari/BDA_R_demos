#' ---
#' title: "Bayesian data analysis demo 6.1"
#' author: "Aki Vehtari, Markus Paasiniemi"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     theme: readable
#'     code_download: true
#' ---

#' ## Posterior predictive checking of normal model for light data
#' 

#' ggplot2 and bayesplot are used for plotting
#' tidyr is used for manipulating data frames
#' posterior is used to handle posterior draws
#+ setup, message=FALSE, error=FALSE, warning=FALSE
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(posterior)
library(tidyr)
library(latex2exp)
library(rprojroot)
root<-has_file(".BDA_R_demos_root")$make_fix_file()

#' Data
y <- read.table(root("demos_ch6","light.txt"))$V1
#' Sufficient statistics
n <- length(y)
s <- sd(y)
my <- mean(y)

#' Create 9 random replicate data sets from the posterior
#' predictive density.
#' Each set has same number of virtual observations as the
#' original data set.
yrep <- replicate(9, rt(n, n-1)*sqrt(1+1/n)*s+my)

#' Replace one of the replicates with observed data.
#' If you can spot which one has been replaced, it means
#' that the replicates do not resemble the original data
#' and thus the model has a defect
ind <- sample(9, 1)
yrep_df <- yrep %>%
  as.data.frame() %>%
  replace(ind, y) %>% # replace one of the yreps
  setNames(1:9) %>%   # rename to hide which one is the real data
  pivot_longer(everything()) # use the long format for plotting

ggplot(data = yrep_df) +
  geom_histogram(aes(x = value), fill = 'steelblue',
                 color = 'black', binwidth = 4) +
  facet_wrap(~name, nrow = 3) +
  coord_cartesian(xlim = c(-55, 65)) +
  labs(x = '', y = '') +
  scale_y_continuous(breaks=NULL) +
  theme(strip.background = element_blank())

#' Generate 1000 replicate data sets and compute test statistic. The
#' test statistic here is the smallest observation in the data or in
#' the replicated data.
yrep1000 <- replicate(1000, rt(n, n-1)*sqrt(1+1/n)*s+my) %>%
  as.data.frame()
# the minimum value over 1000 replicates
minvals <- data.frame(x = sapply(yrep1000, min))
#' Plot test statistic for the data and the replicated data sets
title1 <- 'Smallest observation in the replicated
data (hist.) vs in the original data (vertical line)'
ggplot(data = minvals) +
  geom_histogram(aes(x = x), fill = 'steelblue',
                 color = 'black', binwidth = 2) +
  geom_vline(aes(xintercept = min(x)), data = data.frame(x = y),
             color = 'red') +
  coord_cartesian(xlim = c(-50, 20)) +
  labs(x = TeX('Minimum of \\textit{y} and \\textit{y}$^{rep}$'),
       y = '', title = title1) +
  scale_y_continuous(breaks=NULL)

#' Posterior predictive checks provided by bayesplot<br>
#' https://mc-stan.org/bayesplot/reference/PPC-distributions.html
#'

#' For convenience switch to use draws object
rownames(yrep1000) <- paste0("yrep[", 1:66, "]")
yrep_draws <- as_draws_matrix(t(yrep1000))

#' Histogram of y + 8 yrep histograms
ppc_hist(y, yrep_draws[1:8,])

#' Kernel density estimate of y + 100 yrep kernel density estimates
ppc_dens_overlay(y, yrep_draws[1:100,])

#' ECDF of y + 100 yrep ECDFs
ppc_ecdf_overlay(y, yrep_draws[1:100,])

#' Scatterplot of yrep vs y
ppc_scatter(y, yrep_draws[1:4,])+geom_abline()

#' The distribution of a (test) statistic T(yrep) compared to the
#' observed value T(y) computed from the data y. The default test
#' statistic mean is a bad statistic as the model has a parameter for
#' the mean.
ppc_stat(y, yrep_draws)

#' The distribution of a (test) statistic T(yrep) compared to the
#' observed value T(y) computed from the data y. Min and max are often
#' good test statistics for continuous outcomes.
ppc_stat(y, yrep_draws, stat="min")
ppc_stat(y, yrep_draws, stat="max")

#' Show 2 test statistics in one plot
color_scheme_set("brewer-Paired")
ppc_stat_2d(y, yrep_draws, stat=c("min","max"))
color_scheme_set()
