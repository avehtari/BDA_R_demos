# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

library(ggplot2)
library(gridExtra)
library(tidyr)

# Normal approximaton for Bioassay model.

# Bioassay data, (BDA3 page 86)
df1 <- data.frame(
  x = c(-0.86, -0.30, -0.05, 0.73),
  n = c(5, 5, 5, 5),
  y = c(0, 1, 3, 5)
)

# compute the posterior density in grid
#  - usually should be computed in logarithms!
#  - with alternative prior, check that range and spacing of A and B
#    are sensible

A = seq(-4, 8, length.out = 100)
B = seq(-10, 40, length.out = 100)

# make vectors that contain all pairwise combinations of A and B
cA <- rep(A, each = length(B))
cB <- rep(B, length(A))

# a helper function to calculate the log likelihood
logl <- function(df, a, b)
  df['y']*(a + b*df['x']) - df['n']*log1p(exp(a + b*df['x']))

# calculate likelihoods: apply logl function for each observation
# ie. each row of data frame of x, n and y
p <- apply(df1, 1, logl, cA, cB) %>% rowSums() %>% exp()


# sample from the grid (with replacement)
nsamp <- 1000
samp_indices <- sample(length(p), size = nsamp,
                       replace = T, prob = p/sum(p))

samp_A <- cA[samp_indices[1:nsamp]]
samp_B <- cB[samp_indices[1:nsamp]]
# add random jitter, see BDA3 p. 76
samp_A <- samp_A + runif(nsamp, A[1] - A[2], A[2] - A[1])
samp_B <- samp_B + runif(nsamp, B[1] - B[2], B[2] - B[1])

# samples of LD50 conditional beta > 0
bpi <- samp_B > 0
samp_ld50 <- -samp_A[bpi]/samp_B[bpi]

# limits for the plots
xl <- c(-2, 8)
yl <- c(-10, 40)

# create a plot of the posterior density
pos <- ggplot(data = data.frame(cA ,cB, p), aes(x = cA, y = cB)) +
  geom_raster(aes(fill = p, alpha = p), interpolate = T) +
  geom_contour(aes(z = p), colour = 'black', size = 0.2) +
  coord_cartesian(xlim = xl, ylim = yl) +
  labs(x = 'alpha', y = 'beta') +
  scale_fill_gradient(low = 'yellow', high = 'red', guide = F) +
  scale_alpha(range = c(0, 1), guide = F)

# plot of the samples
sam <- ggplot(data = data.frame(samp_A, samp_B)) +
  geom_point(aes(samp_A, samp_B), color = 'blue', size = 0.3) +
  coord_cartesian(xlim = xl, ylim = yl) +
  labs(x = 'alpha', y = 'beta')

# plot of the histogram of LD50
his <- ggplot() +
  geom_histogram(aes(samp_ld50), binwidth = 0.04,
                 fill = 'steelblue', color = 'black') +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  labs(x = 'LD50 = -alpha/beta')

# define the function to be optimized
bioassayfun <- function(w, df) {
  z <- w[1] + w[2]*df$x
  -sum(df$y*(z) - df$n*log1p(exp(z)))
}


# initial guess
w0 <- c(0,0)
optim_res <- optim(w0, bioassayfun, gr = NULL, df1, hessian = T)
w <- optim_res$par
S <- solve(optim_res$hessian)

# multivariate normal probability density function
dmvnorm <- function(x, mu, sig)
  exp(-0.5*(length(x)*log(2*pi) + log(det(sig)) + (x-mu)%*%solve(sig, x-mu)))

# evaluate likelihood at points (cA,cB)
p <- apply(cbind(cA, cB), 1, dmvnorm, w, S)

# sample from the grid (with replacement)
nsamp <- 1000
samp_indices <- sample(length(p), size = nsamp,
                       replace = T, prob = p/sum(p))

samp_A <- cA[samp_indices[1:nsamp]]
samp_B <- cB[samp_indices[1:nsamp]]
# add random jitter, see BDA3 p. 76
samp_A <- samp_A + runif(nsamp, A[1] - A[2], A[2] - A[1])
samp_B <- samp_B + runif(nsamp, B[1] - B[2], B[2] - B[1])

# samples of LD50 conditional beta > 0
# Normal approximation does not take into account that the posterior
# is not symmetric and that there is very low density for negative
# beta values. Based on the draws from the normal approximation
# is is estimated that there is about 6% probability that beta is negative!
bpi <- samp_B > 0
samp_ld50 <- -samp_A[bpi]/samp_B[bpi]


# create a plot of the posterior density
pos_norm <- ggplot(data = data.frame(cA ,cB, p), aes(x = cA, y = cB)) +
  geom_raster(aes(fill = p, alpha = p), interpolate = T) +
  geom_contour(aes(z = p), colour = 'black', size = 0.2) +
  coord_cartesian(xlim = xl, ylim = yl) +
  labs(x = 'alpha', y = 'beta') +
  scale_fill_gradient(low = 'yellow', high = 'red', guide = F) +
  scale_alpha(range = c(0, 1), guide = F)

# plot of the samples
sam_norm <- ggplot(data = data.frame(samp_A, samp_B)) +
  geom_point(aes(samp_A, samp_B), color = 'blue', size = 0.3) +
  coord_cartesian(xlim = xl, ylim = yl) +
  labs(x = 'alpha', y = 'beta')

# plot of the histogram of LD50
his_norm <- ggplot() +
  geom_histogram(aes(samp_ld50), binwidth = 0.04,
                 fill = 'steelblue', color = 'black') +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  labs(x = 'LD50 = -alpha/beta, beta > 0')

# combine the plots
grid.arrange(pos, sam, his, pos_norm, sam_norm, his_norm, ncol = 3)

