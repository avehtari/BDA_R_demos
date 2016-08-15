# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

library(ggplot2)
library(tidyr)

# Binomial example - Testing sequential dependence example
# Testing sequential dependence example (Gelman et al p. 163)
y <- c(1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
Ty <- sum(diff(y) != 0) + 0.0

# sufficient statistics
n <- length(y)
s <- sum(y)

rb <- function(s, n) {
  p <- rbeta(1, s+1, n-s+1)
  yr <- rbinom(n, 1, p)
  sum(diff(yr) != 0) + 0.0
}

Tyr <- data.frame(x = replicate(10000, rb(s, n)))
mean(Tyr<=Ty)

title1 <- 'Binomial example - number of changes?
Pr(T(yrep,theta) <= T(y,theta)|y) = 0.03'
ggplot(data = Tyr) +
  geom_histogram(aes(x = x), fill = 'darkblue',
                 color = 'black', binwidth = 1) +
  geom_vline(aes(xintercept = x), data = data.frame(x = Ty),
             color = 'red') +
  labs(x = '', y = '', title = title1) +
  scale_y_continuous(breaks=NULL)
# vertical line corresponds to the original data,
# histogram to the generated data

