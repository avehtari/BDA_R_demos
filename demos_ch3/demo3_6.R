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

# Illustrate the posterior inference for Bioassay data (BDA3 p. 74-)

# data
x <- c(-0.86, -0.30, -0.05, 0.73)
n <- c(5, 5, 5, 5)
y <- c(0, 1, 3, 5)

# compute the posterior density in grid
#  - usually should be computed in logarithms!
#  - with alternative prior, check that range and spacing of A and B
#    are sensible

A = seq(-4, 8, length.out = 50)
B = seq(-10, 40, length.out = 50)

# make vectors that contain all pairwise combinations of A and B
cA <- rep(A,each = length(B))
cB <- rep(B,length(A))



# make a helper function to calculate log likelihood
logl <- function(a,b,x,n,y) y*(a+b*x)-n*log1p(exp(a+b*x))

# calculate likelihoods:
# apply logl function for each observation
p <- apply(data.frame(x,n,y),1,
           function(d) logl(cA,cB,d[1],d[2],d[3])) %>%
  # sum the log likelihoods of observations
  # and exponentiate to get the joint likelihood
  rowSums() %>% exp()

# sample from the grid (with replacement)
nsamp <- 1000
samp_indices <- sample.int(n = length(p), size = nsamp,
                       replace = T, prob = p/sum(p))

samp_A <- cA[samp_indices[1:nsamp]]
samp_B <- cB[samp_indices[1:nsamp]]
# add random jitter, see BDA3 p. 76
samp_A <- samp_A + runif(nsamp,A[1]-A[2],A[2]-A[1])
samp_B <- samp_B + runif(nsamp,B[1]-B[2],B[2]-B[1])

# samples of LD50 conditional beta > 0
bpi <- samp_B > 0
samp_ld50 <- -samp_A[bpi]/samp_B[bpi]

# limits for the plots
xl <- c(-2,8)
yl <- c(-2,40)

# create a plot of the posterior density
pos <- ggplot(data.frame(alpha=cA,beta=cB,p),aes(alpha, beta, fill = p)) +
  geom_raster(interpolate = T) + guides(fill=FALSE) + coord_cartesian(xlim=xl,ylim=yl)

# plot the samples
sam <- ggplot(data.frame(alpha=samp_A,beta=samp_B)) +
  geom_point(aes(alpha,beta)) + coord_cartesian(xlim=xl,ylim=yl)

# plot the histogram of LD50
his <- ggplot() + coord_cartesian(xlim=c(-0.5,0.5)) + xlab('ld50=-alpha/beta') +
  geom_histogram(aes(samp_ld50),bins=42,fill='steelblue',color='black')

# combine the plots
grid.arrange(pos,sam,his,nrow=3)
