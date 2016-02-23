# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# ggplot2 is used for plotting, gridExtra for showing multiple
# plots side by side and tidyr for manipulating data frames
if(!require(ggplot2)) install.packages('ggplot2'); require(ggplot2)
if(!require(gridExtra)) install.packages('gridExtra'); require(gridExtra)
if(!require(tidyr)) install.packages('tidyr'); require(tidyr)

# fake interesting distribution
x <- seq(-3, 3, length.out = 200)
r <- c(1.1, 1.3, -0.1, -0.7, 0.2, -0.4, 0.06, -1.7,
       1.7, 0.3, 0.7, 1.6, -2.06, -0.74, 0.2, 0.5)
# Estimate the density (named q, to emphesize that it does not need to be
# normalized). Parameter bw=0.5 is used to mimic the outcome of the
# kernelp function in Matlab.
q <- density(r, bw = 0.5, n = 200, from = -3, to = 3)$y

# importance sampling example
g <- dnorm(x)
w <- q/g
rs <- rnorm(100)
rs <- rs[abs(rs) < 3] # remove samples out of the grid
# find nearest point for which the kernel has been evaluated for each sample
rsi <- sapply(rs, function(arg) which.min(abs(arg - x)))
wr <- q[rsi]/dnorm(x[rsi])

df1 <- data.frame(x, q, g) %>% gather(grp, p, -x)

# create a plot of the target and proposal distributions
distr <- ggplot(data = df1) +
  geom_line(aes(x, p, fill = grp, color = grp)) +
  labs(title = 'Target and proposal distributions', x = '', y = '') +
  scale_color_discrete(labels = c('q(theta|y)', 'g(theta)')) +
  theme(legend.position = 'bottom', legend.title = element_blank())

# create a plot of the samples and importance weights
samp <- ggplot() +
  geom_line(aes(x, w, color = '1')) +
  geom_segment(aes(x = x[rsi], xend = x[rsi], y = 0, yend = wr),
               alpha = 0.5, color = 'steelblue') +
  labs(title = 'Samples and importance weights', x = '', y = '') +
  scale_color_manual(values = c('steelblue'), labels = 'q(theta|y)/g(theta)') +
  theme(legend.position = 'bottom', legend.title = element_blank())

# combine the plots
grid.arrange(distr, samp)

