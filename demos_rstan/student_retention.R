#' ---
#' title: "Bayesian data analysis - Student retention"
#' author: "Aki Vehtari"
#' date: "First version 2023-11-22. Last modified `r format(Sys.Date())`."
#' output:
#'   html_document:
#'     fig_caption: yes
#'     toc: TRUE
#'     toc_depth: 2
#'     number_sections: FALSE
#'     toc_float:
#'       smooth_scroll: FALSE
#'     theme: readable
#'     code_download: true
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(cache=FALSE, message=FALSE, error=FALSE, warning=TRUE, out.width='95%')

#' # Introduction
#'
#' This includes the code used to create the models and plots for
#' student retention data used as an example in several BDA course
#' lectures in 2023.
#'
#' #### Load packages
library(brms)
#options(brms.backend = "cmdstanr")
# Using RStan backend for moment matching in LOO
options(brms.backend = "rstan")
library(posterior)
library(tidybayes)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(latex2exp)
library(khroma)
library(RColorBrewer)
theme_set(bayesplot::theme_default(base_family = "sans", base_size=16))

#' # Data
#'
#' Since 2018, there have been 9 assignments in BDA course.  Data are
#' the number of students who submitted each assignment. As the course
#' is not compulsory for most, and it is common that some students
#' register for more courses than they need and during the course may
#' decide to drop-out. Also as the students can submit also in spring,
#' some students may stop submitting in autumn and get back in spring.
#' Although there are also external reasons why students drop out from
#' the course, we are interested in following the retention and hope
#' that the changes in the course would improve the retention. As a
#' starting point we want to analyse whether we can see differences in
#' years.
#' 
# Number of students who returned the first assignment on 2018-2022
nstudents1<-rep(c(176,242,332,301,245),each=9)
# Number of students for each round
nstudents<-c(c(176, 174, 158, 135, 138, 129, 126, 123, 121),
             c(242, 212, 184, 177, 174, 172, 163, 156, 153),
             c(332, 310, 278, 258, 243, 242, 226, 224, 218),
             c(301, 269, 231, 232, 217, 208, 193, 191, 190),
             c(245, 240, 228, 217, 206, 199, 191, 182, 175))
# Proportion of students
propstudents<-c(c(176, 174, 158, 135, 138, 129, 126, 123, 121)/176,
                c(242, 212, 184, 177, 174, 172, 163, 156, 153)/242,
                c(332, 310, 278, 258, 243, 242, 226, 224, 218)/332,
                c(301, 269, 231, 232, 217, 208, 193, 191, 190)/301,
                c(245, 240, 228, 217, 206, 199, 191, 182, 175)/245)
# Year as integers and factors
year <- rep(2018:2022,each=9)
fyear <- factor(year)
# Assignment numbers
assignment <- rep(1:9, 5)
# Tibble
tb <- tibble(assignment, nstudents, nstudents1, propstudents, year, fyear)

# Another tibble for including 2023 first submission numbers
# Number of students who returned the first assignment
nstudents1<-rep(c(176,242,332,301,245,264),each=9)
# Number of students for each round
nstudents<-c(c(176, 174, 158, 135, 138, 129, 126, 123, 121),
                c(242, 212, 184, 177, 174, 172, 163, 156, 153),
                c(332, 310, 278, 258, 243, 242, 226, 224, 218),
                c(301, 269, 231, 232, 217, 208, 193, 191, 190),
                c(245, 240, 228, 217, 206, 199, 191, 182, 175),
                c(264,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA))

propstudents<-c(c(176, 174, 158, 135, 138, 129, 126, 123, 121)/176,
                c(242, 212, 184, 177, 174, 172, 163, 156, 153)/242,
                c(332, 310, 278, 258, 243, 242, 226, 224, 218)/332,
                c(301, 269, 231, 232, 217, 208, 193, 191, 190)/301,
                c(245, 240, 228, 217, 206, 199, 191, 182, 175)/245,
                c(264,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA)/264)
# Year as integers and factors
year <- rep(2018:2023,each=9)
fyear <- factor(year)
# Assignment numbers
assignment <- rep(1:9, 6)
# Tibble
tb2 <- tibble(assignment, nstudents, nstudents1, propstudents, year, fyear)

# Another tibble for including 2023 all submission rounds
# Number of students for each round
nstudents<-c(c(176, 174, 158, 135, 138, 129, 126, 123, 121),
                c(242, 212, 184, 177, 174, 172, 163, 156, 153),
                c(332, 310, 278, 258, 243, 242, 226, 224, 218),
                c(301, 269, 231, 232, 217, 208, 193, 191, 190),
                c(245, 240, 228, 217, 206, 199, 191, 182, 175),
                c(264, 249, 215, 221, 215, 206, 192,  186, 179))

propstudents<-c(c(176, 174, 158, 135, 138, 129, 126, 123, 121)/176,
                c(242, 212, 184, 177, 174, 172, 163, 156, 153)/242,
                c(332, 310, 278, 258, 243, 242, 226, 224, 218)/332,
                c(301, 269, 231, 232, 217, 208, 193, 191, 190)/301,
                c(245, 240, 228, 217, 206, 199, 191, 182, 175)/245,
                c(264, 249, 215, 221, 215, 206, 192, 186, 179)/264)
# Year as integers and factors
year <- rep(2018:2023,each=9)
fyear <- factor(year)
# Assignment numbers
assignment <- rep(1:9, 6)
# Tibble
tb3 <- tibble(assignment, nstudents, nstudents1, propstudents, year, fyear)

#' Plot of submission numbers
ggplot(tb, aes(x=assignment, y=nstudents, group=year, color=as.factor(year-2017))) +
  geom_line(linewidth=1) + geom_point(size=3) +
  scale_x_continuous(breaks=1:9, lim=c(1,10.2)) +
  labs(x="Assignment", y="Number of students submitted", title="Student retention") +
  theme(legend.position="none") +
  scale_color_bright()+
  annotate(geom="text", x=rep(9.9, 6), y=c(121,153,218,192,169,180), label=c(2018:2023), size=5, color = color("bright")(6))

#' Plot of raw proportions *100 for different assignments and years
ggplot(tb, aes(x=assignment, y=propstudents*100, group=year, color=as.factor(year-2017))) +
  geom_line(linewidth=1) + geom_point(size=3) +
  scale_x_continuous(breaks=1:9, lim=c(1,10.2)) +
  ylim(c(55,100)) +
  labs(x="Assignment", y="Proportion submitted %", title="Student retention") +
  theme(legend.position="none") +
  scale_color_bright()+
  annotate(geom="text", x=rep(9.9, 5), y=c(69, 61, 66, 63.5, 72), label=c(2018:2022), size=5, color = color("bright")(5))

#' # Models
#' 
#' Latent hierarchical linear model + binomial observation model
# save_pars is used for later moment matching
#+ results='hide'
fit4 <- brm(nstudents | trials(nstudents1) ~ assignment + (assignment | year), family=binomial(), data=filter(tb, assignment>1), control = list(adapt_delta = 0.95), save_pars=save_pars(all=TRUE), seed=7253)
#' First with plain PSIS-LOO, and we see some warnings
fit4 <- add_criterion(fit4, 'loo', save_psis=TRUE, moment_match=FALSE, overwrite=TRUE)
#+
loo(fit4)
#' PSIS-LOO + moment matching. There is still one khat>0.7, which we
#' could fix with reloo=TRUE, but skip that now.
# overwrite is needed to force LOO recomputation
fit4 <- add_criterion(fit4, 'loo', save_psis=TRUE, moment_match=TRUE, overwrite=TRUE)
loo(fit4)

#' Latent spline + hierarchical linear model with binomial observation model
#+ results='hide'
fit6 <- brm(nstudents | trials(nstudents1) ~ s(assignment, k=4) + (assignment | year), family=binomial(), data=filter(tb, assignment>1), control = list(adapt_delta = 0.95), save_pars=save_pars(all=TRUE), seed=7253)
#' First with plain PSIS-LOO, and we see some warnings
fit6 <- add_criterion(fit6, 'loo', save_psis=TRUE, moment_match=FALSE, overwrite=TRUE)
#+
loo(fit6)
#' PSIS-LOO + moment matching
# overwrite is needed to force LOO recomputation
fit6 <- add_criterion(fit6, 'loo', save_psis=TRUE, moment_match=TRUE, overwrite=TRUE)
loo(fit6)

#' Compare models
loo_compare(loo(fit4), loo(fit6))

#' # Model predictions
#' 
#' Plot intervals for assignment 9 proportion estimates
assign9linpred<-rvar(posterior_linpred(fit6, newdata=filter(tb,assignment==9), trandform=TRUE))
data.frame(year=c("2018","2019","2020","2021","2022"),propstudents=mean(assign9linpred),q05=quantile(assign9linpred,0.05)[1:5],q95=quantile(assign9linpred,0.95)[1:5])|>
ggplot(aes(x=year, y=propstudents*100, ymin=q05*100, ymax=q95*100)) +
  geom_pointrange(color=4) +
  labs(x="Year", y="Proportion submitted %", title="Proportion submitting 9th assign. (90% interval)")

#' Plot distribution of the difference in linear predictor
assign9linpred<-rvar(posterior_linpred(fit6, newdata=filter(tb,assignment==9), transform=TRUE))
(assign9linpred[5]-assign9linpred[1:4]) |>
  as_draws_df() |>
  rename_with(~ c("2018","2019","2020","2021"), starts_with("x")) |>
  pivot_longer(cols=starts_with("2"), names_to="year", values_to="value") |>
  ggplot(aes(y=year, x=value)) +
  stat_halfeye() +
  geom_vline(xintercept=0) +
  ylab("Year")+
  xlab(TeX("2022 retention was better $\\rightarrow$"))

#' Plot model prediction of proportions for 2018-2022
tb |>
  filter(assignment>1)|>
  add_linpred_draws(fit6, transform=TRUE) |>
  ggplot(aes(x = assignment, y = propstudents, group=fyear)) +
  stat_lineribbon(aes(y = .linpred), .width = c(.95), alpha = 1/2, color=brewer.pal(5, "Blues")[[5]])+
  geom_point(data=tb, color=1)+
  scale_fill_brewer() +
  facet_grid(. ~ fyear)+
  theme(legend.position="none")+
  labs(x="Assignment", y = "Proportion of submissions")+
  scale_x_continuous(breaks=1:9, lim=c(1,10.2))

#' Plot model prediction of proportions for 2018-2023 after we have
#' observed only the first assignment submission numbers for 2023
tb2 |>
  filter(assignment>1)|>
  group_by(fyear) |>
  add_linpred_draws(fit6, transform=TRUE, allow_new_levels=TRUE) |>
  ggplot(aes(x = assignment, y = propstudents, group=fyear)) +
  stat_lineribbon(aes(y = .linpred), .width = c(.95), alpha = 1/2, color=brewer.pal(5, "Blues")[[5]])+
  geom_point(data=tb, color=1)+
  scale_fill_brewer() +
  facet_grid(. ~ fyear)+
  theme(legend.position="none")+
  labs(x="Assignment", y = "Proportion of submissions")+
  scale_x_continuous(breaks=1:9, lim=c(1,10.2))

#' Update the posterior with 2023 data
#+ results='hide'
fit6b <- update(fit6, newdata=filter(tb3, assignment>1))

#' Plot model prediction of proportions for 2018-2022 after observing
#' all submission numbers for 2023
tb3 |>
  filter(assignment>1)|>
  add_linpred_draws(fit6b, transform=TRUE) |>
  ggplot(aes(x = assignment, y = propstudents, group=fyear)) +
  stat_lineribbon(aes(y = .linpred), .width = c(.95), alpha = 1/2, color=brewer.pal(5, "Blues")[[5]])+
  geom_point(data=tb3, color=1)+
  scale_fill_brewer() +
  facet_grid(. ~ fyear)+
  theme(legend.position="none")+
  labs(x="Assignment", y = "Proportion of submissions")+
  scale_x_continuous(breaks=1:9, lim=c(1,10.2))

#' The posterior predictive distribution for year 2023
posterior_linpred(fit6, allow_new_levels=TRUE,
                  newdata=filter(tb3,year==2023&assignment==9),
                  transform = TRUE) |>
  quantile(c(0.05,0.95))

#' # PPC
#'
#' PPC density overlays
pp_check(fit4, ndraws=20)+
  xlim(c(50,370))
pp_check(fit6, ndraws=20)+
  xlim(c(50,370))

#' PPC intervals grouped
pp_check(fit4, type = "intervals_grouped", group="year",
         facet_args=list(nrow=1))
pp_check(fit6, type = "intervals_grouped", group="year",
         facet_args=list(nrow=1))

#' PPC ribbon grouped
pp_check(fit4, type = "ribbon_grouped", group="year",
         facet_args=list(nrow=1))
pp_check(fit6, type = "ribbon_grouped", group="year",
         facet_args=list(nrow=1))

#' PPC intervals
pp_check(fit4, type = "intervals")
pp_check(fit6, type = "intervals")

#' PPC LOO intervals. We get LOO warnings, but in this case they don't matter.
pp_check(fit4, type = "loo_intervals")
pp_check(fit6, type = "loo_intervals")

#' PPC LOO-PIT-QQ-uniform
pp_check(fit4, type = "loo_pit_qq", ndraws=4000)
pp_check(fit6, type = "loo_pit_qq", ndraws=4000)

#' PPC LOO-PIT-QQ-normal
pp_check(fit4, type = "loo_pit_qq", ndraws=4000, compare="normal")
pp_check(fit6, type = "loo_pit_qq", ndraws=4000, compare="normal")
