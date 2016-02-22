# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# R versions of demos 3_1-3_4 combined into one demo
# like iPython notebook versions

# ggplot2 is used for plotting
if(!require(ggplot2)) install.packages('ggplot2')
require(ggplot2)
# gridExtra is for showing multiple plots side by side
if(!require(gridExtra)) install.packages('gridExtra')
require(gridExtra)
# tidyr is for manipulating the data frames
if(!require(tidyr)) install.packages('tidyr')
require(tidyr)

# data
y <- c(93,112,122,135,122,150,118,90,124,114)
# sufficient statistics
n <- length(y)
s2 <- var(y)
my <- mean(y)

# Factorize the joint posterior p(mu,sigma2|y) to p(sigma2|y)p(mu|sigma2,y)
# Sample from the joint posterior using this factorization

# helper function to sample from and evaluate
# scaled inverse chi-squared distribution
rsinvchisq <- function(n,nu,s2,...) nu*s2/rchisq(n,nu,...)
dsinvchisq <- function(x,nu,s2){
  exp(log(nu/2)*nu/2 - lgamma(nu/2) + log(s2)/2*nu - log(x)*(nu/2+1) - (nu*s2/2)/x)
}

# sample 1000 random numbers from p(sigma2|y)
sigma2  <- rsinvchisq(1000,n-1,s2)
# sample from p(mu|sigma2,y)
mu <- my+sqrt(sigma2/n)*rnorm(length(sigma2))
# display sigma instead of sigma2
sigma <- sqrt(sigma2)

# sample from the predictive distribution p(ynew|y)
# for each sample of (mu, sigma)
ynew <- rnorm(1000,mu,sigma)

# For mu compute the density in a grid
t1l <- c(90,150)
t2l <- c(10,60)
t1 <- seq(t1l[1], t1l[2], length.out = 1000)
# For sigma compute the density in these points
t2 <- seq(t2l[1],t2l[2], length.out = 1000)

# combine these into data frame
dd <- data.frame(mu = rep(t1,each = length(t2)),
                 sigma = rep(t2,length(t1)))

# evaluate the joint density in a grid
# note that the following is not normalized, but for plotting
# contours it does not matter
dd <- within(dd,z <- dsinvchisq(sigma^2,n-1,s2)*2*sigma*dnorm(mu,my,sigma/sqrt(n)))

# For ynew compute the density in these points
nl = c(50, 185)
xynew <- seq(nl[1], nl[2], length.out = 1000)

# compute the exact marginal density of mu
# multiplication by 1./sqrt(s2/n) is due to the transformation of
# variable z=(x-mean(y))/sqrt(s2/n), see Gelman et al p. 24
pm <- dt((t1-mean(y))/sqrt(s2/n),n-1)/sqrt(s2/n)

# estimate the marginal density using samples
# and ad hoc Gaussian kernel approximation
pmk <- density(mu,adjust=2,n=1000,from=t1l[1],to=t1l[2])$y

# compute the exact marginal density of sigma
# the multiplication by 2*t2 is due to the transformation of
# variable z=t2^2, see Gelman et al p. 24
ps <- dsinvchisq(t2^2,n-1,s2)*2*t2

# estimate the marginal density using samples
# and ad hoc Gaussian kernel approximation
psk <- density(sigma,n=1000,from=t2l[1],to=t2l[2])$y

# compute the exact predictive density
# multiplication by 1./sqrt(s2/n) is due to the transformation of variable
# see BDA3 p. 21
p_new <- dt((xynew-my)/sqrt(s2*(1+1/n)), n-1) / sqrt(s2*(1+1/n))

## Visualise the joint density and marginal densities of the posterior
## of normal distribution with unknown mean and variance.

# create a plot of the marginal density of mu
margmu <- data.frame(mu=t1,Exact=pm,Empirical=pmk) %>% gather(grp,p,Exact,Empirical) %>%
  ggplot(aes(mu,p,color=grp)) +  geom_line() +
  ggtitle('Marginal of mu') + ylab('') + xlab('') +
  scale_y_continuous(breaks=NULL) + xlim(t1l) +
  theme(legend.background=element_blank(),
        legend.position=c(0.75,0.8),
        legend.title = element_blank())

# create a plot of the marginal density of sigma
margsig <- data.frame(sigma=t2,Empirical=psk,Exact=ps) %>% gather(grp,p,Exact,Empirical) %>%
  ggplot(aes(sigma,p,color=grp)) + geom_line() + coord_flip() + ggtitle('Marginal of sigma') +
  ylab('') + xlab('') + scale_y_continuous(breaks=NULL) + xlim(t2l) +
  theme(legend.background=element_blank(), legend.position=c(0.75,0.8), legend.title = element_blank())

# create a plot of the joint density
joint1 <- ggplot() + coord_cartesian(xlim=t1l,ylim=t2l) + xlab('') + ylab('') +
  geom_point(data=data.frame(mu,sigma),aes(mu,sigma,z=NA),size=0.1) +
  stat_contour(data=dd, aes(mu, sigma, z=z),breaks=seq(1e-5,max(dd$z),length.out=6)) +
  scale_y_continuous(labels=NULL) + scale_x_continuous(labels=NULL) + ggtitle('Joint posterior')

# blank plot for combining the plots
bp <- grid.rect(gp=gpar(col='white'))

# combine the plots
grid.arrange(joint1, margsig, margmu, bp,ncol=2, nrow=2)

## Visualise factored sampling and the corresponding
## marginal and conditional densities.

# data frame for the conditional of mu and marginal of sigma
dm <- data.frame(mu = t1, marg = rep(sigma[1],length(t1)),
                 cond = sigma[1] + dnorm(t1,my,sqrt(sigma2[1]/n))*100) %>%
  gather(grp,p,marg,cond)

# legend labels for the following plot
joint2labs <- c('Exact contour plot','Sample from joint post.',
               'Cond. distribution of mu','Sample from the marg. of sigma')

# create another plot of the joint posterior,
# currently uses an awkward hack to construct the legend
joint2 <- ggplot() + coord_cartesian(xlim=t1l,ylim=t2l)  + ggtitle('Joint posterior') +
  stat_contour(data=dd, aes(x=mu, y=sigma, z=z,color='1'),breaks=seq(1e-5,max(dd$z),length.out=6)) +
  geom_point(data=data.frame(x=mu[1],y=sigma[1]),aes(x,y,color='2')) +
  geom_line(data=dm,aes(x=mu,y=p,color=grp),linetype='dashed') +
  guides(color=guide_legend(nrow=2,title=NULL,override.aes=list(shape=c(NA,16,NA,NA),linetype=c(1,0,2,2)))) +
  theme(legend.position = c(0.5,0.85),legend.background=element_blank()) +
  scale_color_manual(values=c('blue', 'black', 'darkred', 'black'), labels=joint2labs)

# create another plot of the marginal density of sigma
margsig2 <- data.frame(sigma = t2,Empirical = psk,Exact = ps) %>% gather(grp,p,Exact,Empirical) %>%
  ggplot(aes(sigma,p,color=grp)) + geom_line() + coord_flip() + ggtitle('Marginal of sigma') + ylab('') +
  theme(legend.background=element_blank(), legend.position=c(0.75,0.8), legend.title = element_blank())

# combine the plots
grid.arrange(joint2,margsig2,ncol=2)

## Visualise the marginal distribution of mu as a mixture of normals.

# calculate conditional pdfs for each sample
condpdfs <- sapply(t1,function(x) dnorm(x,my,sqrt(sigma2/n)))

# create a plot of some of them
condmu <- data.frame(mu=t1,t(condpdfs[1:25,])) %>% gather(obs,p,X1:X25) %>%
  ggplot() + geom_line(aes(mu,p,group=obs),linetype='dashed') + ylab('') +
  scale_y_continuous(breaks=NULL) + ggtitle('conditional distribution of mu for first 25 samples')

# labels for the following plot
mulabs <- c('avg of sampled conds','exact marginal of mu')

# create a plot of their mean
meanmu <- data.frame(mu = t1, avg = colMeans(condpdfs), exact = pm) %>% gather(grp,p,avg,exact) %>%
  ggplot() + geom_line(aes(mu,p,size=grp,color=grp)) + scale_size_manual(values=c(2,0.8),labels=mulabs) +
  scale_color_manual(values=c('orange','black'),labels=mulabs) + scale_y_continuous(breaks=NULL) + ylab('') +
  theme(legend.position = c(0.8,0.8),legend.background=element_blank(),legend.title = element_blank())

# combine the plots
grid.arrange(condmu,meanmu,ncol=1)

## Visualise sampling from the posterior predictive distribution.

# calculate each predictive pdf with given mu and sigma
ynewdists <- sapply(xynew,function(x) dnorm(x, mu, sigma))

# legend labels for the plots
joint3labs <- c('exact contour plot','samples')
pred1labs <- c('Sample from the predictive distribution','Predictive distribution given the posterior sample')
pred2labs <- c('Samples from the predictive distribution', 'Exact predictive distribution')

# create yet another plot of the joint posterior with a draw
# from the posterior predictive distribution, highlight the first sample
# currently uses an awkward hack to construct the legend
joint3 <- ggplot() + coord_cartesian(xlim=t1l,ylim=t2l)  + xlab('mu') + ylab('sigma') +
  geom_point(data=data.frame(x=mu,y=sigma),aes(x,y,color='2'),size=0.1) +
  stat_contour(data=dd, aes(x=mu, y=sigma, z=z, color='1'),breaks=seq(1e-5,max(dd$z),length.out=6)) +
  geom_point(data=data.frame(x=mu[1],y=sigma[1]),aes(x,y),color='red') +
  guides(color=guide_legend(nrow=2,override.aes=list(linetype=c(1,0),shape=c(NA,16)))) +
  scale_color_manual(values=c('blue','black'),labels=joint3labs) +
  theme(legend.position='bottom', legend.title = element_blank())

# create a plot of the predicitive distribution and the respective sample,
# currently uses an awkward hack to construct the legend
pred1 <- ggplot(data=data.frame(ytilde=xynew,p=ynewdists[1,])) + ylab('') +
  geom_line(aes(ytilde,p,color='2')) + scale_y_continuous(breaks=NULL) +
  geom_point(data=data.frame(ytilde=ynew[1],p=0.02*max(ynewdists)),aes(ytilde,p,color='1')) +
  guides(color=guide_legend(nrow=2,override.aes=list(linetype=c(0,1),shape=c(16,NA),labels = pred1labs))) +
  scale_color_manual(values=c('red','blue'),labels=pred1labs) + xlim(nl) +
  theme(legend.position='bottom', legend.title = element_blank())

# create a plot for all ynews
pred2 <- ggplot(data=data.frame(ytilde=xynew,p=p_new)) + ylab('') +
  geom_line(aes(ytilde,p,color='2')) + scale_y_continuous(breaks=NULL) +
  geom_point(data=data.frame(ytilde=ynew,p=0.02*max(ynewdists)),aes(ytilde,p,color='1'),alpha=0.1) +
  guides(color=guide_legend(nrow=2,override.aes=list(linetype=c(0,1),shape=c(16,NA),labels = pred1labs))) +
  scale_color_manual(values=c('darkred','blue'),labels=pred2labs) +
  theme(legend.position='bottom', legend.title = element_blank()) +
  coord_cartesian(xlim=nl)

# combine the plots
grid.arrange(joint3,pred1,bp,pred2,nrow=2)
