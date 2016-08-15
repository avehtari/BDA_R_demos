# Bayesian data analysis
# Aki Vehtari <Aki.Vehtari@aalto.fi>
# Markus Paasiniemi <Markus.Paasiniemi@aalto.fi>

# Stan demos

library(tidyr)
library(rstan) # version >= 2.10
library(loo)
library(ggplot2)
library(gridExtra)

# note that it is often more convenient to define
# stan models in a separate .stan-file instead of
# a character vector in r

# # # #
# Bernoulli model
code_bern <- '
data {
  int<lower=0> N;
  int<lower=0,upper=1> y[N];
}
parameters {
  real<lower=0,upper=1> theta;
}
model {
  theta ~ beta(1,1);
  for (n in 1:N)
    y[n] ~ bernoulli(theta);
}'

d_bern <- list(N = 10, y = c(0, 1, 0, 0, 1, 1, 1, 0, 1, 0))

fit_bern <- stan(model_code = code_bern, data = d_bern)
print(fit_bern)
samples_bern <- extract(fit_bern, permuted = T)
qplot(x = samples_bern$theta, geom = 'histogram', bins = 50)
# or simply plot(fit_bern)


# # # #
# vectorized Bernoulli model
code_bern <- '
data {
  int<lower=0> N;
  int<lower=0,upper=1> y[N];
}
parameters {
  real<lower=0,upper=1> theta;
}
model {
  theta ~ beta(1,1);
  y ~ bernoulli(theta);
}'

fit_bern <- stan(model_code = code_bern, data = d_bern)
print(fit_bern)
samples_bern <- extract(fit_bern, permuted = T)
qplot(x = samples_bern$theta, geom = 'histogram', bins = 50)

# # # #
# Binomial model
code_bin <- '
data {
  int<lower=0> N;
  int<lower=0> y;
}
parameters {
  real<lower=0,upper=1> theta;
}
model {
  theta ~ beta(1,1);
  y ~ binomial(N,theta);
}'

d_bin <- list(N = 10, y = 8)
fit_bin <- stan(model_code = code_bin, data = d_bin)
print(fit_bin)
samples_bin <- extract(fit_bin, permuted = T)
qplot(x = samples_bin1$theta, geom = 'histogram', bins = 50)

# re-running Binomial model with new data
d_bin <- list(N = 10, y = 10)
fit_bin <- stan(model_code = code_bin, data = d_bin)
print(fit_bin)
samples_bin <- extract(fit_bin, permuted = T)
qplot(x = samples_bin$theta, geom = 'histogram', bins = 50)

# # # #
# comparison of two groups with Binomial
code_bin <- '
data {
  int<lower=0> N1;
  int<lower=0> y1;
  int<lower=0> N2;
  int<lower=0> y2;
}
parameters {
  real<lower=0,upper=1> theta1;
  real<lower=0,upper=1> theta2;
}
model {
  theta1 ~ beta(1,1);
  theta2 ~ beta(1,1);
  y1 ~ binomial(N1,theta1);
  y2 ~ binomial(N2,theta2);
}
generated quantities {
  real oddsratio;
  oddsratio = (theta2/(1-theta2))/(theta1/(1-theta1));
}'

d_bin <- list(N1 = 674, y1 = 39, N2 = 680, y2 = 22)
fit_bin <- stan(model_code = code_bin, data = d_bin)
print(fit_bin)
samples_bin <- extract(fit_bin, permuted = T)
qplot(x = samples_bin$oddsratio, geom = 'histogram', bins = 50)

# # # #
# Gaussian linear model
code_lin <- '
data {
  int<lower=0> N; // number of data points
  vector[N] x; //
  vector[N] y; //
  real xpred; // input location for prediction
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] mu;
  mu = alpha + beta*x;
}
model {
  y ~ normal(mu, sigma);
}
generated quantities {
  real ypred;
  vector[N] log_lik;
  ypred = normal_rng(alpha + beta*xpred, sigma);
  for(i in 1:N)
    log_lik[i] = normal_lpdf(y[i] | mu[i], sigma);
}'

d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_lin <- list(x = rep(d_kilpis$year, each = 3), xpred = 2016,
              y = c(t(d_kilpis[,2:4])), N = 3*nrow(d_kilpis))

# compile and fit the model
fit_lin <- stan(model_code = code_lin, data = d_lin)
samples_lin <- extract(fit_lin, permuted = T)

# plot
mu <- apply(samples_lin$mu, 2, quantile, c(0.05, 0.5, 0.95)) %>%
  t() %>% data.frame(x = d_lin$x, .)  %>% gather(pct, y, -x)
ypred <- samples_lin$ypred
beta <- samples_lin$beta
mean(beta>0) # probability that beta > 0

p1 <- ggplot() +
  geom_point(aes(x, y), data = data.frame(d_lin), size = 0.5) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(y = 'Summer temp. @Kilpisjärvi') +
  guides(linetype = F)
p2 <- qplot(x = beta, geom = 'histogram', bins = 50)
p3 <- qplot(x = ypred, geom = 'histogram', bins = 50) +
  labs(x = 'y-prediction for x=2016')
grid.arrange(p1, p2, p3, nrow = 1)

# psis-loo
log_lik <- extract_log_lik(fit_lin)
loo1 <- loo(log_lik)

# # # #
# Gaussian linear model with adjustable priors
code_lin <- '
data {
  int<lower=0> N; // number of data points
  vector[N] x; //
  vector[N] y; //
  real pmualpha; // prior mean for alpha
  real psalpha;  // prior std for alpha
  real pmubeta;  // prior mean for beta
  real psbeta;   // prior std for beta
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] mu;
  mu = alpha + beta*x;
}
model {
  alpha ~ normal(pmualpha, psalpha);
  beta ~ normal(pmubeta, psbeta);
  y ~ normal(mu, sigma);
}'

d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_lin <-list(N = 3*nrow(d_kilpis),
  x = rep(d_kilpis$year, each = 3),
  y = c(t(d_kilpis[,2:4])),
  pmualpha = mean(unlist(d_kilpis[,2:4])), # centered
  psalpha = (14-4)/6, # avg temp between 4-14
  pmubeta = 0, # a priori incr. and decr. as likely
  psbeta = (.1--.1)/6) # avg temp prob does does not
                       # incr. more than a degree per 10 years

# compile and fit the model
fit_lin <- stan(model_code = code_lin, data = d_lin)
samples_lin <- extract(fit_lin, permuted = T)

# plot
mu <- apply(samples_lin$mu, 2, quantile, c(0.05, 0.5, 0.95)) %>%
  t() %>% data.frame(x = d_lin$x, .)  %>% gather(pct, y, -x)
sigma <- samples_lin$sigma
beta <- samples_lin$beta
mean(beta>0) # probability that beta > 0

p1 <- ggplot() +
  geom_point(aes(x, y), data = data.frame(d_lin), size = 0.5) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(y = 'Summer temp. @Kilpisjärvi') +
  guides(linetype = F)
p2 <- qplot(x = beta, geom = 'histogram', bins = 50)
p3 <- qplot(x = sigma, geom = 'histogram', bins = 50)
grid.arrange(p1, p2, p3, nrow = 1)

# # # #
# Gaussian linear model with standardized data
# this is alternative to above
code_lin <- '
data {
  int<lower=0> N; // number of data points
  vector[N] x; //
  vector[N] y; //
}
transformed data {
  vector[N] x_std;
  vector[N] y_std;
  x_std = (x - mean(x)) / sd(x);
  y_std = (y - mean(y)) / sd(y);
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma_std;
}
transformed parameters {
  vector[N] mu_std;
  mu_std = alpha + beta*x_std;
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  y_std ~ normal(mu_std, sigma_std);
}
generated quantities {
  vector[N] mu;
  real<lower=0> sigma;
  mu = mean(y) + mu_std*sd(y);
  sigma = sigma_std*sd(y);
}'

d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_lin <-list(N = 3*nrow(d_kilpis),
  x = rep(d_kilpis$year, each = 3),
  y = c(t(d_kilpis[,2:4])))

# compile and fit the model
fit_lin <- stan(model_code = code_lin, data = d_lin)
samples_lin <- extract(fit_lin, permuted = T)

# plot
mu <- apply(samples_lin$mu, 2, quantile, c(0.05, 0.5, 0.95)) %>%
  t() %>% data.frame(x = d_lin$x, .)  %>% gather(pct, y, -x)
sigma <- samples_lin$sigma
beta <- samples_lin$beta
mean(beta>0) # probability that beta > 0

p1 <- ggplot() +
  geom_point(aes(x, y), data = data.frame(d_lin), size = 0.5) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(y = 'Summer temp. @Kilpisjärvi') +
  guides(linetype = F)
p2 <- qplot(x = beta, geom = 'histogram', bins = 50)
p3 <- qplot(x = sigma, geom = 'histogram', bins = 50)
grid.arrange(p1, p2, p3, nrow = 1)

# # # #
# Gaussian linear student-t model
code_lin <- '
data {
  int<lower=0> N; // number of data points
  vector[N] x; //
  vector[N] y; //
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
  real<lower=1, upper=80> nu;
}
transformed parameters {
  vector[N] mu;
  mu = alpha + beta*x;
}
model {
  nu ~ gamma(2, 0.1); // Juárez and Steel(2010)
  y ~ student_t(nu, mu, sigma);
}'

d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_lin <-list(N = 3*nrow(d_kilpis),
  x = rep(d_kilpis$year, each = 3),
  y = c(t(d_kilpis[,2:4])))

# compile and fit the model
fit_lin <- stan(model_code = code_lin, data = d_lin)
samples_lin <- extract(fit_lin, permuted = T)

# plot
mu <- apply(samples_lin$mu, 2, quantile, c(0.05, 0.5, 0.95)) %>%
  t() %>% data.frame(x = d_lin$x, .)  %>% gather(pct, y, -x)
sigma <- samples_lin$sigma
beta <- samples_lin$beta
mean(beta>0) # probability that beta > 0
nu <- samples_lin$nu

p1 <- ggplot() +
  geom_point(aes(x, y), data = data.frame(d_lin), size = 0.5) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(y = 'Summer temp. @Kilpisjärvi') +
  guides(linetype = F)
p2 <- qplot(x = beta, geom = 'histogram', bins = 50)
p3 <- qplot(x = sigma, geom = 'histogram', bins = 50)
p4 <- qplot(x = nu, geom = 'histogram', bins = 50)
grid.arrange(p1, p2, p3, p4)

# # # #
# comparison of k groups with common variance (ANOVA)
code_grp <- '
data {
  int<lower=0> N; // number of data points
  int<lower=0> K; // number of groups
  int<lower=1,upper=K> x[N]; // group indicator
  vector[N] y; //
}
parameters {
  vector[K] mu;        // group means
  real<lower=0> sigma; // common std
}
model {
  for (n in 1:N)
    y[n] ~ normal(mu[x[n]], sigma);
}'

d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_grp <-list(N = 3*nrow(d_kilpis), K = 3,
             x = rep(1:3, nrow(d_kilpis)),
             y = c(t(d_kilpis[,2:4])))

fit_grp <- stan(model_code = code_grp, data = d_grp)
samples_grp <- extract(fit_grp, permuted = T)

mu <- data.frame(samples_grp$mu) %>% setNames(6:8) %>% gather(month, temp)
qplot(month, temp, data = mu, geom = 'boxplot')

# probabilities that june is hotter than july, june is hotter than august
# and july is hotter than august:
combn(unique(mu$month), 2, function(a, df)
  mean(subset(df, month == a[1])$temp > subset(df, month == a[2])$temp), df = mu)

# # # #
# comparison of k groups with unequal variances
code_grp <- '
data {
  int<lower=0> N; // number of data points
  int<lower=0> K; // number of groups
  int<lower=1,upper=K> x[N]; // group indicator
  vector[N] y; //
}
parameters {
  vector[K] mu;        // group means
  vector<lower=0>[K] sigma; // group stds
}
model {
  for (n in 1:N)
    y[n] ~ normal(mu[x[n]], sigma[x[n]]);
}'

d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_grp <-list(N = 3*nrow(d_kilpis), K = 3,
             x = rep(1:3, nrow(d_kilpis)),
             y = c(t(d_kilpis[,2:4])))

fit_grp <- stan(model_code = code_grp, data = d_grp)
samples_grp <- extract(fit_grp, permuted = T)

mu <- data.frame(samples_grp$mu) %>% setNames(6:8) %>% gather(month, temp)
qplot(month, temp, data = mu, geom = 'boxplot')

# probabilities that june is hotter than july, june is hotter than august
# and july is hotter than august:
combn(unique(mu$month), 2, function(a, df)
  mean(subset(df, month == a[1])$temp > subset(df, month == a[2])$temp), df = mu)

# # # #
# Hierarchical prior for means in comparison of k groups
# results do not differ much from the previous, because there is only
# few groups and quite much data per group, but this works as an example anyway
code_grp <- '
data {
  int<lower=0> N; // number of data points
  int<lower=0> K; // number of groups
  int<lower=1,upper=K> x[N]; // group indicator
  vector[N] y; //
}
parameters {
  real mu0;             // prior mean
  real<lower=0> sigma0; // prior std
  vector[K] mu;         // group means
  real<lower=0> sigma;  // common std
}
model {
  mu0 ~ normal(10, 10);     // weakly informative prior
  sigma0 ~ cauchy(0, 4);    // weakly informative prior
  mu ~ normal(mu0, sigma0); // population prior with unknown parameters
  sigma ~ cauchy(0, 4);     // weakly informative prior
  for (n in 1:N)
    y[n] ~ normal(mu[x[n]], sigma);
}'

d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_grp <-list(N = 3*nrow(d_kilpis), K = 3,
             x = rep(1:3, nrow(d_kilpis)),
             y = c(t(d_kilpis[,2:4])))

fit_grp <- stan(model_code = code_grp, data = d_grp)
samples_grp <- extract(fit_grp, permuted = T)

mu0 <- samples_grp$mu
sd(mu0)
mu <- data.frame(samples_grp$mu) %>% setNames(6:8) %>% gather(month, temp)
qplot(month, temp, data = mu, geom = 'boxplot')

# probabilities that june is hotter than july, june is hotter than august
# and july is hotter than august:
combn(unique(mu$month), 2, function(a, df)
  mean(subset(df, month == a[1])$temp > subset(df, month == a[2])$temp), df = mu)

# # # #
# Hierarchical prior for means in comparison of k groups
# results do not differ much from the previous, because there is only
# few groups and quite much data per group, but this works as an example anyway
code_grp <- '
data {
  int<lower=0> N; // number of data points
  int<lower=0> K; // number of groups
  int<lower=1,upper=K> x[N]; // group indicator
  vector[N] y; //
}
parameters {
  real mu0;                 // prior mean
  real<lower=0> musigma0;   // prior std
  vector[K] mu;             // group means
  real lsigma0;             // prior mean
  real<lower=0> lsigma0s;   // prior std
  vector<lower=0>[K] sigma; // group stds
}
model {
  mu0 ~ normal(10, 10);       // weakly informative prior
  musigma0 ~ cauchy(0,10);    // weakly informative prior
  mu ~ normal(mu0, musigma0); // population prior with unknown parameters
  lsigma0 ~ normal(0,1);      // weakly informative prior
  lsigma0s ~ normal(0,1);     // weakly informative prior
  sigma ~ cauchy(lsigma0, lsigma0s); // population prior with unknown parameters
  for (n in 1:N)
    y[n] ~ normal(mu[x[n]], sigma[x[n]]);
}'

d_kilpis <- read.delim('kilpisjarvi-summer-temp.csv', sep = ';')
d_grp <-list(N = 3*nrow(d_kilpis), K = 3,
             x = rep(1:3, nrow(d_kilpis)),
             y = c(t(d_kilpis[,2:4])))

fit_grp <- stan(model_code = code_grp, data = d_grp)
samples_grp <- extract(fit_grp, permuted = T)

mu0 <- samples_grp$mu
sd(mu0)
mu <- data.frame(samples_grp$mu) %>% setNames(6:8) %>% gather(month, temp)
qplot(month, temp, data = mu, geom = 'boxplot')

# probabilities that june is hotter than july, june is hotter than august
# and july is hotter than august:
combn(unique(mu$month), 2, function(a, df)
  mean(subset(df, month == a[1])$temp > subset(df, month == a[2])$temp), df = mu)
