%   Author: Vehtari Aki <Aki.Vehtari@aalto.fi>
%   Last modified: 2015-12-09 16:34:21 EET

% When running in brute.aalto.fi
addpath ~ave/matlab/MatlabProcessManager
addpath ~ave/matlab/MatlabStan

%% Bernoulli model
bernoulli_code={
  'data {'
  '  int<lower=0> N;'
  '  int<lower=0,upper=1> y[N];'
  '}'
  'parameters {'
  '  real<lower=0,upper=1> theta;'
  '}'
  'model {'
  '  theta ~ beta(1,1);'
  '  for (n in 1:N)'
  '    y[n] ~ bernoulli(theta);'
  '}'
  };
dat=struct('N',10,'y',[0 1 0 0 1 1 1 0 1 0]);
fit = stan('model_code',bernoulli_code,'data',dat,'sample_file','bernoulli','file_overwrite',true,'verbose',true);
fit.block()
print(fit);
samples = fit.extract('permuted',true);
hist(samples.theta,50)

%% Vectorized Bernoulli model
bernoulli_code={
  'data {'
  '  int<lower=0> N;'
  '  int<lower=0,upper=1> y[N];'
  '}'
  'parameters {'
  '  real<lower=0,upper=1> theta;'
  '}'
  'model {'
  '  theta ~ beta(1,1);'
  '  y ~ bernoulli(theta);'
  '}'
  };
dat=struct('N',10,'y',[1 1 1 0 1 1 1 0 1 1]);
fit = stan('model_code',bernoulli_code,'data',dat,'sample_file','bernoulli','file_overwrite',true,'verbose',true);
fit.block()

%% Binomial model
binomial_code={
  'data {'
  '  int<lower=0> N;'
  '  int<lower=0> y;'
  '}'
  'parameters {'
  '  real<lower=0,upper=1> theta;'
  '}'
  'model {'
  '  theta ~ beta(1,1);'
  '  y ~ binomial(N,theta);'
  '}'
              };
dat=struct('N',10,'y',8);
fit = stan('model_code',binomial_code,'data',dat,'sample_file','binomial','file_overwrite',true,'verbose',true);
fit.block()
samples = fit.extract('permuted',true);
hist(samples.theta,50)

%% re-running Binomial model with new data
dat=struct('N',10,'y',10);
fit = stan('fit',fit,'data',dat,'sample_file','binomial','file_overwrite',true,'verbose',false);
fit.block()
samples = fit.extract('permuted',true);
hist(samples.theta,50)

%% Comparison of two groups with Binomial 
binomial_code={
  'data {'
  '  int<lower=0> N1;'
  '  int<lower=0> y1;'
  '  int<lower=0> N2;'
  '  int<lower=0> y2;'
  '}'
  'parameters {'
  '  real<lower=0,upper=1> theta1;'
  '  real<lower=0,upper=1> theta2;'
  '}'
  'model {'
  '  theta1 ~ beta(1,1);'
  '  theta2 ~ beta(1,1);'
  '  y1 ~ binomial(N1,theta1);'
  '  y2 ~ binomial(N2,theta2);'
  '}'
  'generated quantities {'
  '  real oddsratio;'
  '  oddsratio <- (theta2/(1-theta2))/(theta1/(1-theta1));'
  '}'
              };
dat=struct('N1',674,'y1',39,'N2',680,'y2',22);
fit = stan('model_code',binomial_code,'data',dat,'sample_file','binomial','file_overwrite',true,'verbose',true);
fit.block()
samples = fit.extract('permuted',true);
hist(samples.oddsratio,50)

%% Gaussian linear model
linear_code = {
   'data {'
   '    int<lower=0> N; // number of data points '
   '    vector[N] x; // '
   '    vector[N] y; // '
   '    real xpred; // input location for prediction'
   '}'
   'parameters {'
   '    real alpha; '
   '    real beta; '
   '    real<lower=0> sigma;'
   '}'
   'transformed parameters {'
   '    vector[N] mu;'
   '    mu <- alpha + beta*x;'
   '}'
   'model {'
   '    y ~ normal(mu, sigma);'
   '}'
   'generated quantities {'
   '    real ypred;'
   '    vector[N] log_lik;'
   '    ypred <- normal_rng(alpha + beta*xpred, sigma);'
   '    for (n in 1:N)'
   '        log_lik[n] <- normal_log(y[n], alpha + beta*x[n], sigma);'
   '}'
};
% Data for Stan
d = dataset('File','kilpisjarvi-summer-temp.csv','Delimiter',';','ReadVarNames',true);
x = repmat(double(d(:,1)),1,3)'; x = x(:);
y = double(d(:,2:4))'; y = y(:);
N = numel(x);
xpred = 2016;
dat = struct('N',N, 'x',x,'y',y,'xpred',xpred);
% Compile and fit the model
fit = stan('model_code',linear_code,'data',dat,'sample_file','kilpis','file_overwrite',true,'verbose',true);
fit.block()

% Plot
samples=fit.extract('permuted',true);
subplot(1,3,1)
mu = samples.mu;
ypred = samples.ypred;
plot(x,prctile(mu,[50]),'r-',x,prctile(mu,[5
95]),'r--',x,y,'b.')
xlabel('Year')
ylabel('Summer temperature @ Kilpisjarvi');
subplot(1,3,2)
beta = samples.beta;
hist(beta,50)
xlabel('beta')
%probability that beta>0
mean(beta>0)
subplot(1,3,3)
hist(samples.ypred,50)
xlabel(sprintf('y-prediction for x=%d', xpred))

% psis-loo
[loo,~,khat] = psisloo(samples.log_lik);

%% Gaussian linear model with adjustable priors
linear_code = {
  'data {'
  '    int<lower=0> N; // number of data points '
  '    vector[N] x; // '
  '    vector[N] y; // '
  '    real pmualpha; // prior mean for alpha'
  '    real psalpha;  // prior std for alpha'
  '    real pmubeta;  // prior mean for beta'
  '    real psbeta;   // prior std for beta'
  '}'
  'parameters {'
  '    real alpha; '
  '    real beta; '
  '    real<lower=0> sigma;'
  '}'
  'transformed parameters {'
  '    vector[N] mu;'
  '    mu <- alpha + beta*x;'
  '}'
  'model {'
  '    alpha ~ normal(pmualpha,psalpha);'
  '    beta ~ normal(pmubeta,psbeta);'
  '    y ~ normal(mu, sigma);'
  '}'
};
% Data for Stan
d=dataset('File','kilpisjarvi-summer-temp.csv','Delimiter',';','ReadVarNames',true);
x=repmat(double(d(:,1)),1,3)';x=x(:);
y=double(d(:,2:4))';y=y(:);
N=numel(x);
dat = struct('N',N, ...
             'x',x, ...
             'y',y, ...
             'pmualpha',mean(y), ... % centered
             'psalpha',(14-4)/6, ... % avg temp between 4-14
             'pmubeta',0,...        % a priori increase and decrese as likely
             'psbeta',(.1--.1)/6 ... % avg temp probably does not increase more than 1 degree per 10 years
             );
% Compile and fit the model
fit = stan('model_code',linear_code,'data',dat,'sample_file','kilpis','file_overwrite',true,'verbose',true);
fit.block()

% Plot
samples=fit.extract('permuted',true);
subplot(1,3,1)
mu = samples.mu;
plot(x,prctile(mu,[50]),'r-',x,prctile(mu,[5 95]),'r--',x,y,'b.')
xlabel('Year')
ylabel('Summer temperature @ Kilpisjarvi');
subplot(1,3,2)
beta = samples.beta;
hist(beta,50)
xlabel('beta')
%probability that beta>0
mean(beta>0)
subplot(1,3,3)
sigma = samples.sigma;
hist(sigma,50)
xlabel('sigma')

%% Gaussian linear model with standardized data
% This is alternative to above
linear_code = {
   'data {'
   '    int<lower=0> N; // number of data points '
   '    vector[N] x; // '
   '    vector[N] y; // '
   '}'
  'transformed data {'
  '  vector[N] x_std;'
  '  vector[N] y_std;'
  '  x_std <- (x - mean(x)) / sd(x);'
  '  y_std <- (y - mean(y)) / sd(y);'
  '}'
  'parameters {'
  '    real alpha; '
  '    real beta; '
  '    real<lower=0> sigma_std;'
  '}'
  'transformed parameters {'
  '    vector[N] mu_std;'
  '    mu_std <- alpha + beta*x_std;'
  '}'
  'model {'
  '  alpha ~ normal(0,1);'
  '  beta ~ normal(0,1);'
  '  y_std ~ normal(mu_std, sigma_std);'
  '}'
  'generated quantities {'
  '    vector[N] mu;'
  '    real<lower=0> sigma;'
  '    mu <- mean(y) + mu_std*sd(y);'
  '    sigma <- sigma_std*sd(y);'
  '}'
};
% Data for Stan
d=dataset('File','kilpisjarvi-summer-temp.csv','Delimiter',';','ReadVarNames',true);
x=repmat(double(d(:,1)),1,3)';x=x(:);
y=double(d(:,2:4))';y=y(:);
N=numel(x);
dat = struct('N',N,...
             'x',x,...
             'y',y);
% Compile and fit the model
fit = stan('model_code',linear_code,'data',dat,'sample_file','kilpis','file_overwrite',true,'verbose',true);
fit.block()

% Plot
samples=fit.extract('permuted',true);
subplot(1,3,1)
mu = samples.mu;
plot(x,prctile(mu,[50]),'r-',x,prctile(mu,[5 95]),'r--',x,y,'b.')
xlabel('Year')
ylabel('Summer temperature @ Kilpisjarvi');
subplot(1,3,2)
beta = samples.beta;
hist(beta,50)
xlabel('beta')
%probability that beta>0
mean(beta>0)
subplot(1,3,3)
sigma = samples.sigma;
hist(sigma,50)
xlabel('sigma')

%% Gaussian linear student-t model
linear_code = {
   'data {'
   '    int<lower=0> N; // number of data points '
   '    vector[N] x; // '
   '    vector[N] y; // '
   '}'
   'parameters {'
   '    real alpha; '
   '    real beta; '
   '    real<lower=0> sigma;'
   '    real<lower=1,upper=80> nu;'
   '}'
   'transformed parameters {'
   '    vector[N] mu;'
   '    mu <- alpha + beta*x;'
   '}'
   'model {'
   '    nu ~ gamma(2,0.1); // Juárez and Steel (2010)'
   '    y ~ student_t(nu, mu, sigma);'
   '}'
};
% Data for Stan
d=dataset('File','kilpisjarvi-summer-temp.csv','Delimiter',';','ReadVarNames',true);
x=repmat(double(d(:,1)),1,3)';x=x(:);
y=double(d(:,2:4))';y=y(:);
N=numel(x);
dat = struct('N',N,...
             'x',x,...
             'y',y);
% Compile and fit the model
fit = stan('model_code',linear_code,'data',dat,'sample_file','kilpis','file_overwrite',true,'verbose',true);
fit.block()

% Plot
samples=fit.extract('permuted',true);
subplot(2,2,1)
mu = samples.mu;
plot(x,prctile(mu,[50]),'r-',x,prctile(mu,[5 95]),'r--',x,y,'b.')
xlabel('Year')
ylabel('Summer temperature @ Kilpisjarvi');
subplot(2,2,2)
beta = samples.beta;
hist(beta,50)
xlabel('beta')
mean(beta>0)
subplot(2,2,3)
sigma = samples.sigma;
hist(sigma,50)
xlabel('sigma')
subplot(2,2,4)
nu = samples.nu;
hist(nu,50)
xlabel('nu')

%% Comparison of k groups with common variance (ANOVA)
group_code = {
  'data {'
  '    int<lower=0> N; // number of data points '
  '    int<lower=0> K; // number of groups '
  '    int<lower=1,upper=K> x[N]; // group indicator '
  '    vector[N] y; // '
  '}'
  'parameters {'
  '    vector[K] mu;        // group means '
  '    real<lower=0> sigma; // common std '
  '}'
  'model {'
  '    for (n in 1:N)'
  '      y[n] ~ normal(mu[x[n]], sigma);'
  '}'
             };
% Data for Stan
d=dataset('File','kilpisjarvi-summer-temp.csv','Delimiter',';','ReadVarNames',true);
% Is there difference between different summer months?
x=repmat([1:3]',size(d,1),1); % summer months are numbered from 1 to 3
y=double(d(:,2:4))';y=y(:);
N=numel(x);
dat = struct('N',N,...
             'K',3,... % 3 groups
             'x',x,... % group indicators
             'y',y);   % observations
% Compile and fit the model
fit = stan('model_code',group_code,'data',dat,'sample_file','kilpis','file_overwrite',true,'verbose',true);
fit.block()

% Plot
samples=fit.extract('permuted',true);
mu = samples.mu;
boxplot(mu)
% matrix of probabilities that one mu is larger than other
for k1=1:3
  for k2=(k1+1):3
    ps(k1,k2)=mean(mu(:,k1)>mu(:,k2));
    ps(k2,k1)=1-ps(k1,k2);
  end
end
ps

%% Comparison of k groups with unequal variances
group_code = {
  'data {'
  '    int<lower=0> N; // number of data points '
  '    int<lower=0> K; // number of groups '
  '    int<lower=1,upper=K> x[N]; // group indicator '
  '    vector[N] y; // '
  '}'
  'parameters {'
  '    vector[K] mu;    // group means '
  '    vector<lower=0>[K] sigma; // group stds '
  '}'
  'model {'
  '    for (n in 1:N)'
  '      y[n] ~ normal(mu[x[n]], sigma[x[n]]);'
  '}'
             };
% Data for Stan
d=dataset('File','kilpisjarvi-summer-temp.csv','Delimiter',';','ReadVarNames',true);
% Is there difference between different summer months?
x=repmat([1:3]',size(d,1),1); % summer months are numbered from 1 to 3
y=double(d(:,2:4))';y=y(:);
N=numel(x);
dat = struct('N',N,...
             'K',3,... % 3 groups
             'x',x,... % group indicators
             'y',y);   % observations
% Compile and fit the model
fit = stan('model_code',group_code,'data',dat,'sample_file','kilpis','file_overwrite',true,'verbose',true);
fit.block()

% Plot
samples=fit.extract('permuted',true);
mu = samples.mu;
boxplot(mu)
% matrix of probabilities that one mu is larger than other
for k1=1:3
  for k2=(k1+1):3
    ps(k1,k2)=mean(mu(:,k1)>mu(:,k2));
    ps(k2,k1)=1-ps(k1,k2);
  end
end
ps

%% Hierarchical prior for means in comparison of k groups
% results do not differ much from the previous, because there is only
% few groups and quite much data per group, but this works as an example anyway
hier_code = {
  'data {'
  '    int<lower=0> N; // number of data points '
  '    int<lower=0> K; // number of groups '
  '    int<lower=1,upper=K> x[N]; // group indicator '
  '    vector[N] y; // '
  '}'
  'parameters {'
  '    real mu0;             // prior mean '
  '    real<lower=0> sigma0; // prior std '
  '    vector[K] mu;         // group means '
  '    real<lower=0> sigma;  // common std '
  '}'
  'model {'
  '    mu0 ~ normal(10,10);      // weakly informative prior '
  '    sigma0 ~ cauchy(0,4);     // weakly informative prior '
  '    mu ~ normal(mu0, sigma0); // population prior with unknown parameters'
  '    sigma ~ cauchy(0,4);      // weakly informative prior '
  '    for (n in 1:N)'
  '      y[n] ~ normal(mu[x[n]], sigma);'
  '}'
             };
% Data for Stan
d=dataset('File','kilpisjarvi-summer-temp.csv','Delimiter',';','ReadVarNames',true);
% Is there difference between different summer months?
x=repmat([1:3]',size(d,1),1); % summer months are numbered from 1 to 3
y=double(d(:,2:4))';y=y(:);
N=numel(x);
dat = struct('N',N,...
             'K',3,... % 3 groups
             'x',x,... % group indicators
             'y',y);   % observations
% Compile and fit the model
fit = stan('model_code',hier_code,'data',dat,'sample_file','kilpis','file_overwrite',true,'verbose',true);
fit.block()

% Plot
samples=fit.extract('permuted',true);
mu0 = samples.mu0;
std(mu0)
mu = samples.mu;
boxplot(mu)
% matrix of probabilities that one mu is larger than other
for k1=1:4
  for k2=(k1+1):4
    ps(k1,k2)=mean(mu(:,k1)>mu(:,k2));
    ps(k2,k1)=1-ps(k1,k2);
  end
end
ps

%% Hierarchical prior for means and variances in comparison of k groups
% results do not differ much from the previous, because there is only
% few groups and quite much data per group, but this works as an example anyway
hier_code = {
  'data {'
  '    int<lower=0> N; // number of data points '
  '    int<lower=0> K; // number of groups '
  '    int<lower=1,upper=K> x[N]; // group indicator '
  '    vector[N] y; // '
  '}'
  'parameters {'
  '    real mu0;               // prior mean '
  '    real<lower=0> musigma0; // prior std '
  '    vector[K] mu;           // group means '
  '    real lsigma0;            '
  '    real<lower=0> lsigma0s;  '
  '    vector<lower=0>[K] sigma; // group stds '
  '}'
  'model {'
  '    mu0 ~ normal(10,10);         // weakly informative prior '
  '    musigma0 ~ cauchy(0,10);     // weakly informative prior '
  '    mu ~ normal(mu0, musigma0);  // population prior with unknown parameters'
  '    lsigma0 ~ normal(0,1);       // weakly informative prior '
  '    lsigma0s ~ normal(0,1);      // weakly informative prior '
  '    sigma ~ lognormal(lsigma0, lsigma0s); // population prior with unknown parameters'
  '    for (n in 1:N)'
  '      y[n] ~ normal(mu[x[n]], sigma[x[n]]);'
  '}'
             };
% Data for Stan
d=dataset('File','kilpisjarvi-summer-temp.csv','Delimiter',';','ReadVarNames',true);
% Is there difference between different summer months?
x=repmat([1:3]',size(d,1),1); % summer months are numbered from 1 to 3
y=double(d(:,2:4))';y=y(:);
N=numel(x);
dat = struct('N',N,...
             'K',3,... % 3 groups
             'x',x,... % group indicators
             'y',y);   % observations
% Compile and fit the model
fit = stan('model_code',hier_code,'data',dat,'sample_file','kilpis','file_overwrite',true,'verbose',true);
fit.block()

% Plot
samples=fit.extract('permuted',true);
mu0 = samples.mu0;
std(mu0)
mu = samples.mu;
boxplot(mu)
% matrix of probabilities that one mu is larger than other
for k1=1:4
  for k2=(k1+1):4
    ps(k1,k2)=mean(mu(:,k1)>mu(:,k2));
    ps(k2,k1)=1-ps(k1,k2);
  end
end
ps


