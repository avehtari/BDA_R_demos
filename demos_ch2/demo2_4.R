# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# ggplot2 is used for plotting
if(!require(ggplot2)) install.packages("ggplot2")
require(ggplot2)
# gridExtra is for showing multiple plots side by side
if(!require(gridExtra)) install.packages("gridExtra")
require(gridExtra)
# tidyr is for manipulating the data frames
if(!require(tidyr)) install.packages("tidyr")
require(tidyr)

# Calculate the posterior distribution on a discrete grid of points by
# multiplying the likelihood and a non-conjugate prior at each point,
# and normalizing over the points. Simulate samples from the resulting
# non-standard posterior distribution using inverse cdf using the
# discrete grid.

# Posterior with uniform prior (Beta(1,1))
# data (437,543)
a <- 437
b <- 543
x <- seq(0.1,1,0.001)
con <- dbeta(x,a,b)

# Compute the density of non-conjugate prior in discrete points, i.e. in a grid
# this non-conjugate prior is the same as in figure 2.4 in the book
pp <- rep(1,length(x))
pi1 <- which(x==0.385)
pi2 <- which(x==0.485)
pi3 <- which(x==0.585)
pm <- 11
pp[pi1:pi2] <- seq(1, pm, length.out = length(pi1:pi2))
pp[pi3:pi2] <- seq(1, pm, length.out = length(pi3:pi2))
nc_p <- pp/sum(pp)

# compute the un-normalized non-conjugate posterior in a grid
po <- dbeta(x,a,b)*pp

# normalize the posterior
nc_po <- po/sum(po)

# add variables to data frame
d <- data.frame(x,con,nc_p,nc_po)

# plot posterior with uniform prior, non-conjugate
# prior and the corresponding non-conjugate posterior
gather(d,group,p,con,nc_p,nc_po) %>%
  within(group<-factor(group, labels=c("Posterior with uniform prior",
                                       "Non-conjugate prior",
                                       "Non-conjugate posterior"))) %>%
  ggplot(aes(x=x,y=p)) + geom_line() +
  facet_wrap(~group,ncol=1,scales="free_y") +
  coord_cartesian(xlim=c(0.35,0.6)) +
  ylab("") + scale_y_continuous(breaks=NULL)


# next demonstrate inverse cdf sampling

# compute the cumulative density in a grid
pc <- cumsum(nc_po)
qc <- c(0,pc)

# runif(k) returns k uniform random numbers from interval [0,1]
# set.seed(seed) is used to set seed foor the randon number generator
set.seed(2601)
r <- runif(10000)

# invcdf returns the values of x at points for cdf(r0)<=x
# ie. invcdf(0.5) returns the x-value of the median
invcdf <- function(r0){x[sum(qc<r0)]}

# sapply(x,fun) applies function fun to each element of x
# ie. elements of s are now draws from the distribution
s <- sapply(r,invcdf)

# construct plots: p1 is the posterior, p2 is the cdf of the posterior
# and p3 is the histogram of posterior samples (drawn using inv. cdf)
p1 <- qplot(x,nc_po,geom="line",main="Non-conjugate posterior")
p2 <- qplot(x,pc,geom="line",main="Posterior-cdf")
p3 <- qplot(s,geom="histogram",ylab="",bins=30,main="Histogram of posterior samples")
l <- list(p1,p2,p3)

#write a helper function that adds options to figures
opts <- function(arg){arg + xlab("") +  coord_cartesian(xlim=c(0.35,0.6)) +
    scale_y_continuous(breaks=NULL) + ylab("")}

#apply helper function to each figure
ll <- lapply(l,opts)

# do.call(fun,list) evaluates fun with the elements of list as arguments
# ie. do.call(grid.arrange,list(p1,p2,p3)) => grid.arrange(p1,p2,p3)
# grid.arrange() plots its arguments in one figure
do.call(grid.arrange,ll)
