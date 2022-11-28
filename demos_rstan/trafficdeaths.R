#' ---
#' title: "Bayesian data analysis - traffic deaths in Finland"
#' author: "Aki Vehtari"
#' date: "First version 2017-09-28. Last modified `r format(Sys.Date())`."
#' output:
#'   html_document:
#'     fig_caption: yes
#'     toc: TRUE
#'     toc_depth: 2
#'     number_sections: TRUE
#'     toc_float:
#'       smooth_scroll: FALSE
#'     theme: readable
#'     code_download: true
#' ---

#' # Setup  {.unnumbered}

#+ setup, include=FALSE
knitr::opts_chunk$set(cache=FALSE, message=FALSE, error=FALSE, warning=TRUE, out.width='95%')


#' **Load packages**
#+  comment=NA
library(ggplot2)
library(tidyr)
library(gridExtra)
library(rstanarm)
library(brms)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size = 16))
library(patchwork)
library(loo)
library(rprojroot)
root<-has_file(".BDA_R_demos_root")$make_fix_file()


#' # Introduction
#' 
#' This notebook demonstrates time series analysis for traffic deaths per
#' year in Finland. Currently when the the number of traffic deaths
#' during previous year are reported, the press release claims that the
#' the traffic safety in Finland has improved or worsened depending
#' whether the number is smaller or larger than the year before. Time
#' series analysis can be used to separate random fluctuation from the
#' slowly changing traffic safety.
#' 
#' # Data
#' 
#' Read the data (there would data for earlier years, too, but this is
#' sufficient for the demonstration)
#+ 
# file preview shows a header row
deaths <- read.csv(root("demos_rstan", "trafficdeaths.csv"), header = TRUE)
head(deaths)


#' First plot just the data.
#+
deaths |>
  ggplot(aes(x=year, y=deaths)) +
  geom_point() +
  labs(y = 'Traffic deaths', x= "Year") +
  guides(linetype = "none")

#' # Poisson regression model
#' 
#' The number of deaths is count data, so we use Poisson observation
#' model.  We first fit log-linear model for the Poisson intensity, which
#' corresponds to assuming constant proportional change in the rate.

#+ 
fit_lin <- stan_glm(deaths ~ year, data=deaths, family=poisson,
 	            refresh=1000, iter=1000, chains=4, seed=583829, refresh=0)


#' ESS's and Rhat's are ok (see, e.g., [RStan
#' workflow](http://mc-stan.org/users/documentation/case-studies/rstan_workflow.html)). Let's
#' look at the posterior predictive distribution (median and 5% and 95%
#' intervals).
#+ 
x_predict <- seq(1993,2030)
N_predict <- length(x_predict)
y_predict_lin <- posterior_predict(fit_lin, newdata=data.frame(year=x_predict))
mu <- apply(t(y_predict_lin), 1, quantile, c(0.05, 0.5, 0.95)) %>%
  t() %>% data.frame(x = x_predict, .) %>% gather(pct, y, -x)
pfit <- ggplot() +
  geom_point(aes(year, deaths), data = deaths, size = 1) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(x = 'Year', y = 'Traffic deaths') +
  guides(linetype = F)
(pfit)


#' Next we fit a non-linear spline model with `stan_gamm4`
#+ 
fit_gam <- stan_gamm4(deaths ~ year + s(year), data=deaths,
                      family=poisson, adapt_delta=0.999, 
                      refresh=1000, iter=2000, chain=4, seed=583829, refresh=0)

#' ESS is clearly smaller than for the linear model, but Rhat's are ok.
#' 
#' Let's look at the posterior predictive distribution.
#+ 
x_predict=seq(1993,2030)
N_predict=length(x_predict)
y_predict_gam <- posterior_predict(fit_gam, newdata=data.frame(year=x_predict))
mu <- apply(t(y_predict_gam), 1, quantile, c(0.05, 0.5, 0.95)) %>%
  t() %>% data.frame(x = x_predict, .) %>% gather(pct, y, -x)
pfit <- ggplot() +
  geom_point(aes(year, deaths), data = deaths, size = 1) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(x = 'Year', y = 'Traffic deaths') +
  guides(linetype = F)
(pfit)


#' The predictive median is clearly nonlinear. The predictive mean for
#' future years stays at the same level as the most recent observations,
#' but uncertainty increases quickly.
#' 
#' Finally we fit Gaussian process centered on linear model. We use
#' brms for this:
#+  comment=NA
fit_gp <- brm(deaths ~ year + gp(year), data=deaths,
                      family=poisson, adapt_delta=0.95, 
              refresh=1000, iter=2000, chain=4, seed=583829, refresh=0,
              backend='cmdstanr')

x_predict=seq(1993,2030)
N_predict=length(x_predict)
y_predict_gp <- posterior_predict(fit_gp, newdata=data.frame(year=x_predict))
mu <- apply(t(y_predict_gp), 1, quantile, c(0.05, 0.5, 0.95)) %>%
  t() %>% data.frame(x = x_predict, .) %>% gather(pct, y, -x)
pfit <- ggplot() +
  geom_point(aes(year, deaths), data = deaths, size = 1) +
  geom_line(aes(x, y, linetype = pct), data = mu, color = 'red') +
  scale_linetype_manual(values = c(2,1,2)) +
  labs(x = 'Year', y = 'Traffic deaths') +
  guides(linetype = F)
(pfit)

#' Finally we compare models using PSIS-LOO predictive performance estimates.

#+ 
(loo_lin<-loo(fit_lin, save_psis=TRUE))
(loo_gam<-loo(fit_gam, save_psis=TRUE))
(loo_gp<-loo(fit_gp, save_psis=TRUE))

loo_compare(loo_lin, loo_gam, loo_gp)

#' There are no practical differences in predictive performance, which is
#' partially due to small number of observations. Based on the posterior
#' predictive distributions there are clear differences in the future
#' predictions.

#' We can also look at the calibration of leave-one-out predictive
#' distributions
pi_lin <- ppc_loo_intervals(deaths$deaths,
                    y_predict_lin[,1:29],
                    psis_object=loo_lin$psis_object)+
  labs(title='PPC-LOO linear model')

pi_gp <- ppc_loo_intervals(deaths$deaths,
                    y_predict_gp[,1:29],
                    psis_object=loo_gp$psis_object)+
  labs(title='PPC-LOO GP model')

pi_lin/pi_gp

#' There is a small difference in favor of GP model.

#' <br />
#' 
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2017-2020, Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2017-2020, Aki Vehtari, licensed under CC-BY-NC 4.0.
#' 
#' # Original Computing Environment {.unnumbered}

#+ 
sessionInfo()


#' <br />
