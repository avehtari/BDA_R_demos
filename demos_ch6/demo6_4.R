#' ---
#' title: "Bayesian data analysis demo 6.4"
#' author: "Aki Vehtari, Markus Paasiniemi"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     theme: readable
#'     code_download: true
#' ---

#' ## Marginal posterior predictive checking
#' 
#' Tail area probabilities of marginal predictive distributions,
#' aka probability integral transformation (PIT).
#' 
#' Normal model for light speed data.
#' 

#' ggplot2 is used for plotting, tidyr for manipulating data frames
#+ setup, message=FALSE, error=FALSE, warning=FALSE
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(tidyr)
library(latex2exp)
library(rprojroot)
root<-has_file(".BDA_R_demos_root")$make_fix_file()

#' Data
y <- read.table(root("demos_ch6","light.txt"))$V1
#' Sufficient statistics
n <- length(y)
s <- sd(y)
my <- mean(y)

#' Tail area probabilities of marginal predictive distributions,
#' aka probability integral transformation (PIT)
Ty <- data.frame(x = pt((y - my)/(sqrt(1+1/n)*s), n-1))

#* Plot histogram of PIT values. Ideally histogram should be close to uniform.
title1 <- 'Light speed example
distribution of predictive distribution tail-values'

ggplot(data = Ty) +
  geom_histogram(aes(x = x), fill = 'steelblue',
                 color = 'black', binwidth = 0.05) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = TeX(r"($\textit{p}(\textit{y}^{\textrm{rep}}_{\textit{i}} < \textit{y_i} | \textit{y})$)"),
       y = '', title = title1) +
  scale_y_continuous(breaks=NULL)

#* Plot ECDF of PIT values. Ideally ECDF should be close to diagonal line
ggplot(data=data.frame(x=Ty), aes(x)) +
  stat_ecdf(geom = "step", color=4) +
  xlim(c(0,1))+
  labs(x="Observed PIT values", y="ECDF")+
  annotate(geom="segment",x=0,y=0,xend=1,yend=1)

#' Repeat the PIT checking after removing two "outliers"
y <- y[y>0]
#' Sufficient statistics
n <- length(y)
s <- sd(y)
my <- mean(y)

#' Tail area probabilities of marginal predictive distributions,
#' aka probability integral transformation (PIT)
Ty <- data.frame(x = pt((y - my)/(sqrt(1+1/n)*s), n-1))

#' Plot histogram of PIT values. Ideally histogram should be close to uniform.
title1 <- 'Light speed example
distribution of predictive distribution tail-values'
ggplot(data = Ty) +
  geom_histogram(aes(x = x), fill = 'steelblue',
                 color = 'black', binwidth = 0.05) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = TeX(r"($\textit{p}(\textit{y}^{\textrm{rep}}_{\textit{i}} < \textit{y_i} | \textit{y})$)"),
       y = '', title = title1) +
  scale_y_continuous(breaks=NULL)

#* Plot ECDF of PIT values. Ideally ECDF should be close to diagonal line
ggplot(data=data.frame(x=Ty), aes(x)) +
  stat_ecdf(geom = "step", color=4) +
  xlim(c(0,1))+
  labs(x="Observed PIT values", y="ECDF")+
  annotate(geom="segment",x=0,y=0,xend=1,yend=1)
