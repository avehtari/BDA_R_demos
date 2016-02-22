# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# Illustrate the effect of prior.
# Comparison of posterior distributions with different
# parameter values for the beta prior distribution.

#ggplot2 is used for plotting
if(!require(ggplot2)) install.packages('ggplot2')
require(ggplot2)
#tidyr is used for manipulating data frames
if(!require(tidyr)) install.packages('tidyr')
require(tidyr)

# Prior, posterior and data:
# observed data, 437 girls and 543 boys
a <- 437
b <- 543

# Evaluate densities at evenly spaced points between 0.375 and 0.525
x <- seq(0.375,0.525,0.001)

# Posterior with Beta(1,1), ie. uniform prior
d <- data.frame(x=x,pu=dbeta(x,a+1,b+1))

# 3 different choices for priors, Beta(0.485*2,(1-0.485)*2),
# Beta(0.485*20,(1-0.485)*20) and Beta(0.485*200,(1-0.485)*200)
apr <- 0.485*c(2,20,200)
bpr <- (1-0.485)*c(2,20,200)

# Posteriors
apo <- a + apr
bpo <- b + bpr

figtitles <- paste('alpha/(alpha+beta)=0.485, alpha+beta=',c(2,20,200),sep='')
dfs <- list()
for(i in 1:3){
  # dataframes for each choice of prior for plotting
  dfs[[i]] <- within(d,{figtitle <- figtitles[i]
                        pr <- dbeta(x,apr[i],bpr[i])
                        po <- dbeta(x,apo[i],bpo[i])})
}

# %>% is used for piping
# merge dataframes
Reduce(rbind,dfs) %>%
  # collapse variables into key-value pairs and change their names
  gather(group,p,pu,pr,po) %>%
  within(group<-factor(group, labels=c('Posterior',
                                       'Prior',
                                       'Post with unif prior'))) %>%
  # plot distributions
  ggplot(aes(x=x,y=p,color=group)) + facet_wrap(~figtitle,ncol=1) +
  geom_line(show.legend=T) + ylab('') + xlab('') +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position='bottom',legend.title=element_blank()) +
  geom_vline(xintercept = 0.485,linetype='dotted')

