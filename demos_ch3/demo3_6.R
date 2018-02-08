#' ---
#' title: "Bayesian data analysis demo 3.6"
#' author: "Aki Vehtari, Markus Paasiniemi"
#' date: "`r format(Sys.Date())`"
#' ---

#' ## Illustrate the posterior inference for Bioassay data (BDA3 p. 74-)
#' 

#' ggplot2, and gridExtra are used for plotting, tidyr for manipulating data frames
#+ setup, message=FALSE, error=FALSE, warning=FALSE
library(ggplot2)
library(gridExtra)
library(tidyr)

#' Data
df1 <- data.frame(
  x = c(-0.86, -0.30, -0.05, 0.73),
  n = c(5, 5, 5, 5),
  y = c(0, 1, 3, 5)
)

#' Compute the posterior density in grid
#' 
#'  - usually should be computed in logarithms!
#'  - with alternative prior, check that range and spacing of A and B are sensible
A = seq(-4, 8, length.out = 50)
B = seq(-10, 40, length.out = 50)

#' Make vectors that contain all pairwise combinations of A and B
cA <- rep(A, each = length(B))
cB <- rep(B, length(A))

#' Make a helper function to calculate the log likelihood
#' given a dataframe with x, y, and n and evaluation
#' points a and b. For the likelihood see BDA p. 75
logl <- function(df, a, b)
  df['y']*(a + b*df['x']) - df['n']*log1p(exp(a + b*df['x']))

#' Calculate likelihoods: apply logl function for each observation
#' ie. each row of data frame of x, n and y
p <- apply(df1, 1, logl, cA, cB) %>%
  # sum the log likelihoods of observations
  # and exponentiate to get the joint likelihood
  rowSums() %>% exp()

#' Sample from the grid (with replacement)
nsamp <- 1000
samp_indices <- sample(length(p), size = nsamp,
                       replace = T, prob = p/sum(p))
samp_A <- cA[samp_indices[1:nsamp]]
samp_B <- cB[samp_indices[1:nsamp]]
#' Add random jitter, see BDA3 p. 76
samp_A <- samp_A + runif(nsamp, (A[1] - A[2])/2, (A[2] - A[1])/2)
samp_B <- samp_B + runif(nsamp, (B[1] - B[2])/2, (B[2] - B[1])/2)

#' Sample LD50 conditional beta > 0
bpi <- samp_B > 0
samp_ld50 <- -samp_A[bpi]/samp_B[bpi]

#' Create a plot of the posterior density
# limits for the plots
xl <- c(-2, 8)
yl <- c(-2, 40)
pos <- ggplot(data = data.frame(cA ,cB, p), aes(cA, cB)) +
  geom_raster(aes(fill = p, alpha = p), interpolate = T) +
  geom_contour(aes(z = p), colour = 'black', size = 0.2) +
  coord_cartesian(xlim = xl, ylim = yl) +
  labs(x = 'alpha', y = 'beta') +
  scale_fill_gradient(low = 'yellow', high = 'red', guide = F) +
  scale_alpha(range = c(0, 1), guide = F)

#' Plot of the samples
sam <- ggplot(data = data.frame(samp_A, samp_B)) +
  geom_point(aes(samp_A, samp_B), color = 'blue') +
  coord_cartesian(xlim = xl, ylim = yl) +
  labs(x = 'alpha', y = 'beta')

#' Plot of the histogram of LD50
his <- ggplot() +
  geom_histogram(aes(samp_ld50), binwidth = 0.02,
                 fill = 'steelblue', color = 'black') +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  labs(x = 'LD50 = -alpha/beta')

#' Combine the plots
grid.arrange(pos, sam, his, nrow=3)

