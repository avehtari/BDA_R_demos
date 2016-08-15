# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

library(ggplot2)
library(tidyr)

# Posterior predictive checking
# Light speed example with poorly chosen test statistic

# data
data_path <- 'light.txt'
y <- read.table(data_path)$V1
# sufficient statistics
n <- length(y)
s <- sd(y)
my <- mean(y)

# second example of replications
sampt <- replicate(1000, rt(n, n-1)*sqrt(1+1/n)*s+my) %>%
  as.data.frame()

sampt_vars <- data.frame(x = sapply(sampt, var))

title1 <- 'Light speed example with poorly chosen test statistic
Pr(T(yrep,theta) <= T(y,theta)|y)=0.42'
ggplot(data = sampt_vars) +
  geom_histogram(aes(x = x), fill = 'darkblue',
                 color = 'black', binwidth = 6) +
  geom_vline(aes(xintercept = x), data = data.frame(x = s^2),
             color = 'red') +
  labs(x = '', y = '', title = title1) +
  scale_y_continuous(breaks=NULL)
# vertical line corresponds to the original data,
# histogram to the generated data

