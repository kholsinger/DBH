library(rstan)
library(loo)

rm(list=ls())

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
loo_correlated <- loo(log_lik_gi_correlated)
loo_coupled <- loo(log_lik_gi_coupled)
loo_uncoupled <- loo(log_lik_gi_uncoupled)
print(compare(loo_correlated, loo_coupled, loo_uncoupled), digits=3)

cat("\n\n")
cat("DBH component\n")
loo_correlated <- loo(log_lik_dbh_correlated)
loo_coupled <- loo(log_lik_dbh_coupled)
loo_uncoupled <- loo(log_lik_dbh_uncoupled)
print(compare(loo_correlated, loo_coupled, loo_uncoupled), digits=3)
sink()
