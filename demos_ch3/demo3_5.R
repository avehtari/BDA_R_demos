# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

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
data_path <- 'path/light.txt'
data_path <- 'light.txt'
y <- read.table(data_path)$V1
# sufficient statistics
n <- length(y)
s2 <- var(y)
my <- mean(y)

# filtered data
y_pos <- y[y > 0]
# sufficient statistics
n_pos <- length(y_pos)
s2_pos <- var(y_pos)
my_pos <- mean(y_pos)

# for mu, compute the density in these points
tl1 <- c(18, 34)
t1 <- seq(tl1[1], tl1[2], length.out = 1000)

# compute the exact marginal density for mu
# multiplication by 1./sqrt(s2/n) is due to the transformation of variable
# z=(x-mean(y))/sqrt(s2/n), see BDA3 p. 21
pm_mu <- dt((t1 - my) / sqrt(s2/n), n-1) / sqrt(s2/n)

# compute the exact marginal density for mu for the filtered data
pm_mu_pos = dt((t1 - my_pos) / sqrt(s2_pos/n_pos), n_pos-1) / sqrt(s2_pos/n_pos)

# create a histogram of the measurements
p1 <- ggplot() + ggtitle('Newcomb\'s measurements') +
  geom_histogram(aes(y),bins=42,fill='steelblue',color='black') +
  coord_cartesian(xlim=c(-40,40))

# legend labels for the following plot
labs <- c('posterior of mu', 'posterior of mu given y > 0', '\'true value\'')

# crate a plot of the normal model
p2 <- data.frame(t1,pm_mu,pm_mu_pos) %>%
  gather(grp,p,-t1) %>%
  ggplot() + geom_line(aes(x=t1,y=p,color=grp)) +  ggtitle('Normal model') +
  geom_vline(aes(xintercept = 33, color='intercecpt'),linetype="dashed",show.legend = F) +
  guides(color=guide_legend(override.aes=list(linetype=c(1,1,2),labels = labs))) +
  scale_color_manual(values=c('black','blue','black'),labels=labs) +
  theme(legend.position='bottom', legend.title = element_blank()) +
  coord_cartesian(xlim=c(-40,40))

# combine the plots
grid.arrange(p1,p2,heights = c(2,3))
