# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# ggplot2 is used for plotting
if(!require(ggplot2)) install.packages("ggplot2")
require(ggplot2)
# tidyr is for manipulating the data frames
if(!require(tidyr)) install.packages("tidyr")
require(tidyr)

# Simulate samples from Beta(438,544), draw a histogram with
# quantiles, and do the same for a transformed variable.

# Sample from posterior Beta(438,544)
# Get all samples at once and store them in vector 'th'
a <-438
b <- 544
th <- rbeta(10000,a,b)

x <- seq(0.375,0.525,0.001)

# animation of 200 draws MISSING

# compute 2.5% and 97.5% quantile approximation using samples
thq <- quantile(th,c(0.025,0.975))

# phi = (1-theta)/theta
phi <- (1-th)/th

# 2.5% and 97.5% quantile approximation for phi
phiq <- quantile(phi,c(0.025,0.975))

# merge quantiles into one data frame for plotting
q=gather(data.frame(phi=phiq,th=thq))

# merge the data into one data frame for plotting
dd <- data.frame(phi=phi,th=th)

# Histogram plots with 30 bins for theta and phi
gather(dd) %>%
  ggplot(aes(x=value)) + geom_histogram(bins = 30,fill="darkblue") +
  facet_wrap(~key,ncol=1,scales="free_x")  + ylab("") +
  geom_vline(aes(xintercept = value),data=q,linetype="dotted") +
  scale_y_continuous(breaks=NULL)
