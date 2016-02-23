# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# Rejection sampling example

# ggplot2 is used for plotting, tidyr for manipulating data frames
if(!require(ggplot2)) install.packages('ggplot2'); require(ggplot2)
if(!require(tidyr)) install.packages('tidyr'); require(tidyr)

# fake interesting distribution
x <- seq(-3, 3, length.out = 200)
r <- c(1.1, 1.3, -0.1, -0.7, 0.2, -0.4, 0.06, -1.7,
      1.7, 0.3, 0.7, 1.6, -2.06, -0.74, 0.2, 0.5)
# Estimate the density (named q, to emphesize that it does not need to be
# normalized). Parameter bw=0.5 is used to mimic the outcome of the
# kernelp function in Matlab.
q <- density(r, bw = 0.5, n = 200, from = -3, to = 3)$y

# rejection sampling example
g_mean <- 0
g_sd <- 1.1
g <- dnorm(x, g_mean, g_sd)
# M is computed by discrete approximation
M <- max(q/g)
# prescale
g <- g*M

# illustrate one sample
r1 = -0.8
zi = which.min(abs(x - r1)) # find the closest grid point

r21 = 0.3 * g[zi]
r22 = 0.8 * g[zi]

df1 <- data.frame(x, q, g) %>% gather(grp, p, -x)
# subset with only target distribution
dfq <- subset(df1 , grp == "q")

# labels for the following plot
labs1 <- c('Mg(theta)','q(theta|y)')

# Visualize one accepted and one receted draw:
ggplot() +
  geom_line(data = df1, aes(x, p, fill = grp, color = grp, linetype = grp)) +
  geom_area(data = dfq, aes(x, p), fill = 'lightblue', alpha = 0.3) +
  geom_point(aes(x[zi], r21), col = 'darkgreen', size = 2) +
  geom_point(aes(x[zi], r22), col = 'red', size = 2) +
  geom_segment(aes(x = x[zi], xend = x[zi], y = 0, yend = q[zi])) +
  geom_segment(aes(x = x[zi], xend = x[zi], y = q[zi], yend = g[zi]),
               linetype = 'dashed') +
  scale_y_continuous(breaks = NULL) +
  labs(y = '') +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  scale_linetype_manual(values = c('dashed', 'solid'), labels = labs1) +
  scale_color_discrete(labels = labs1) +
  annotate('text', x = x[zi] + 0.1, y = c(r21, r22),
           label = c('accepted', 'rejected'), hjust = 0)

# get nsamp samples
nsamp <- 200
r1s <- rnorm(nsamp, g_mean, g_sd)
zis <- sapply(r1s, function(r) which.min(abs(x - r)))
r2s <- runif(nsamp) * g[zis]
acc <- ifelse(r2s < q[zis], 'a', 'r')

df2 <- data.frame(r1s, r2s, acc)

#
labs2 <- c('Accepted','Rejected','Mg(theta)','q(theta|y)')

# Visualize 200 draws, only some of which are accepted
ggplot() +
  geom_line(data = df1, aes(x, p, fill = grp, color = grp, linetype = grp)) +
  geom_area(data = dfq, aes(x, p), fill = 'lightblue', alpha = 0.3) +
  geom_point(data = df2, aes(r1s, r2s, color = acc, linetype = acc), size = 2) +
  geom_rug(data = subset(df2, acc== 'a'), aes(x = r1s, r2s),
           col = 'darkgreen', sides = 'b') +
  labs(y = '') +
  scale_y_continuous(breaks = NULL) +
  scale_linetype_manual(values = c(0, 2, 1, 0), labels = labs2) +
  scale_color_manual(values=c('darkgreen','red','black','red'), labels = labs2) +
  guides(color=guide_legend(override.aes=list(
    shape = c(16, 16, NA, NA), linetype = c(0, 0, 2, 1),
    color=c('darkgreen', 'red', 'red', 'black'), labels = labs2))) +
  theme(legend.position = 'bottom', legend.title = element_blank())

# alternative proposal distribution
ga <- rep(0, length(x))
ga[x <= -1.5] <- seq(q[1], max(q[x <= -1.5]), length.out = length(x[x <= -1.5]))
ga[(x > -1.5) & (x <= 0.2)] <- seq(max(q[x <= -1.5]), max(q[(x > -1.5) & (x <= 0.2)]),
                                   length.out = length(x[(x > -1.5) & (x <= 0.2)]))
ga[(x > 0.2) & (x <= 2.3)] <- seq(max(q[(x > -1.5) & (x <= 0.2)]), max(q[x > 2.3]),
                                  length.out = length(x[(x > 0.2) & (x <= 2.3)]))
ga[x > 2.3] <- seq(max(q[x > 2.3]), q[length(q)],
                   length.out = length(x[x > 2.3]))
M = max(q/ga)
ga <- ga*M

df3 <- data.frame(x, q, g = ga) %>% gather(grp, p, -x)

# Visualize alternate proposal distribution
ggplot() +
  geom_line(data = df3, aes(x, p, fill = grp, color = grp, linetype = grp)) +
  geom_area(data = dfq, aes(x, p), fill = 'lightblue', alpha = 0.3) +
  labs(title = 'Alternative proposal distribution', y = '') +
  scale_y_continuous(breaks = NULL) +
  scale_color_discrete(labels = labs1) +
  scale_linetype_manual(values = c('dashed', 'solid'), labels = labs1) +
  theme(legend.position = 'bottom', legend.title = element_blank())

