# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

library(ggplot2)
library(tidyr)

# Posterior predictive checking demo

# data
data_path <- 'light.txt'
y <- read.table(data_path)$V1
# sufficient statistics
n <- length(y)
s <- sd(y)
my <- mean(y)

# Create 10 random replicate data sets from the posterior
# predictive density
# Each set has same number of virtual observations as the
# original data set
sampt <- replicate(10, rt(n, n-1)*sqrt(1+1/n)*s+my) %>%
  as.data.frame()

# Replace one of the replicates with observed data
# If you can spot which one has been replaced, it means
# that the replicates do not resemble the original data
# and thus the model has a defect
ind <- sample(10, 1)
sampt_y <- replace(sampt, ind, y) %>% setNames(1:10) %>% gather()

ggplot(data = sampt_y) +
  geom_histogram(aes(x = value), fill = 'darkblue',
                 color = 'black', binwidth = 5) +
  facet_wrap(~key, nrow = 5) +
  coord_cartesian(xlim = c(-50, 50)) +
  labs(x = '', y = '') +
  scale_y_continuous(breaks=NULL) +
  theme(strip.background = element_blank())


# Generate 1000 replicate data sets
sampt1000 <- replicate(1000, rt(n, n-1)*sqrt(1+1/n)*s+my) %>%
  as.data.frame()
minvals <- data.frame(x = sapply(sampt1000, min))

title1 <- 'Smallest observation in the replicated
data (hist.) vs in the original data (vertical line)'
ggplot(data = minvals) +
  geom_histogram(aes(x = x), fill = 'darkblue',
                 color = 'black', binwidth = 5) +
  geom_vline(aes(xintercept = min(x)), data = data.frame(x = y),
             color = 'red') +
  coord_cartesian(xlim = c(-50, 20)) +
  labs(x = '', y = '', title = title1) +
  scale_y_continuous(breaks=NULL)

