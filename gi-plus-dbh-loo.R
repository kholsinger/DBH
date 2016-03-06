library(rstan)
library(loo)
library(ggplot2)

rm(list=ls())

load("results-gi-plus-dbh-multi-correlated.Rsave")
log_lik_gi_multi_correlated <- extract_log_lik(fit, parameter_name="log_lik_gi")
log_lik_dbh_multi_correlated <- extract_log_lik(fit, parameter_name="log_lik_dbh")

load("results-gi-plus-dbh-correlated.Rsave")
log_lik_gi_correlated <- extract_log_lik(fit, parameter_name="log_lik_gi")
log_lik_dbh_correlated <- extract_log_lik(fit, parameter_name="log_lik_dbh")

load("results-gi-plus-dbh-uncoupled.Rsave")
log_lik_gi_uncoupled <- extract_log_lik(fit, parameter_name="log_lik_gi")
log_lik_dbh_uncoupled <- extract_log_lik(fit, parameter_name="log_lik_dbh")

load("results-gi-plus-dbh.Rsave")
log_lik_gi_coupled <- extract_log_lik(fit, parameter_name="log_lik_gi")
log_lik_dbh_coupled <- extract_log_lik(fit, parameter_name="log_lik_dbh")

sink("model-comparison.txt")
cat("Growth Increment component\n")
loo_multi_correlated <- loo(log_lik_gi_multi_correlated)
loo_correlated <- loo(log_lik_gi_correlated)
loo_coupled <- loo(log_lik_gi_coupled)
loo_uncoupled <- loo(log_lik_gi_uncoupled)
print(compare(loo_multi_correlated, loo_correlated, loo_coupled, loo_uncoupled), digits=3)

n <- dim(loo_correlated$pointwise)[1]
gi <- data.frame(Model=c(rep("Uncoupled", n), rep("Coupled", n), rep("Correlated", n),
                         rep("Multi-correlated", n)),
                 elpd=c(loo_uncoupled$pointwise[,1],
                        loo_coupled$pointwise[,1],
                        loo_correlated$pointwise[,1],
                        loo_multi_correlated$pointwise[,1]),
                 Individual=seq(from=1, to=n, by=1))
p <- ggplot(gi, aes(x=Individual, y=elpd, color=Model)) + geom_point()
print(p)

cat("\n\n")
cat("DBH component\n")
loo_multi_correlated <- loo(log_lik_dbh_multi_correlated)
loo_correlated <- loo(log_lik_dbh_correlated)
loo_coupled <- loo(log_lik_dbh_coupled)
loo_uncoupled <- loo(log_lik_dbh_uncoupled)
print(compare(loo_multi_correlated, loo_correlated, loo_coupled, loo_uncoupled), digits=3)
sink()

n <- dim(loo_correlated$pointwise)[1]
gi <- data.frame(Model=c(rep("Uncoupled", n), rep("Coupled", n), rep("Correlated", n),
                         rep("Multi-correlated", n)),
                 elpd=c(loo_uncoupled$pointwise[,1],
                        loo_coupled$pointwise[,1],
                        loo_correlated$pointwise[,1],
                        loo_multi_correlated$pointwise[,1]),
                 Individual=seq(from=1, to=n, by=1))
p <- ggplot(gi, aes(x=Individual, y=elpd, color=Model)) + geom_point()
print(p)

