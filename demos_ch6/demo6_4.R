# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

library(ggplot2)
library(tidyr)

# Marginal posterior predictive checking
# Light speed example

# data
data_path <- 'light.txt'
y <- read.table(data_path)$V1
# sufficient statistics
n <- length(y)
s <- sd(y)
my <- mean(y)

# tail area probabilities of marginal predictive distributions
Ty <- data.frame(x = pt((y - my)/(sqrt(1+1/n)*s), n-1))

title1 <- 'Light speed example
distribution of marginal posterior p-values'
ggplot(data = Ty) +
  geom_histogram(aes(x = x), fill = 'darkblue',
                 color = 'black', binwidth = 0.05) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = '', y = '', title = title1)

