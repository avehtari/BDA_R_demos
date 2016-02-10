# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

#ggplot2 is used for plotting
if(!require(ggplot2)) install.packages("ggplot2")
require(ggplot2)

# 437 girls and 543 boys have been observed. Calculate and plot the
# posterior distribution of the proportion of girls theta, using
# uniform prior on theta

# posterior is Beta(348,544)
df1 <- data.frame(x=seq(0.375,0.525,0.001))
a <- 438
b <- 544
# dbeta computes the posterior density
df1$p <- dbeta(df1$x,a,b)

# seq creates evenly spaced values from 2.5% quantile to 97.5% quantile (i.e., 95% central interval)
# qbeta computes the value for a given quantile given parameters a and b
df2 <- data.frame(x=seq(qbeta(0.025,a,b),qbeta(0.975,a,b),length.out=100))
# compute the posterior density
df2$p <- dbeta(df2$x,a,b)
df <- bind_rows(df1,df2)

# Plot posterior Beta(438,544)
qplot(x,p,data=df1,geom="line") +
  labs(title="Uniform prior -> Posterior is Beta(438,544)",y="") +
  scale_y_continuous(expand = c(0,0),breaks=NULL) +
  # Add a layer of colorized 95% posterior interval
  geom_area(aes(x,p),data=df2,fill="darkblue",show.legend=T) +
  # Add the proportion of girl babies in general population
  geom_vline(xintercept = 0.485,linetype="dotted")

