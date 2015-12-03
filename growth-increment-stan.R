library(reshape)
library(rstan)

rm(list=ls())

debug <- FALSE

## MCMC settings
##
n.burnin <- 1250
n.iter <- 2500
n.thin <- 1
n.chains <- 4
if (debug) {
  n.chains <- 2
  ## to allow replication of results across runs in JAGS
  ##
}

rstan_options(auto_write = TRUE)
n.cores <- parallel::detectCores()

## select last k rows from x
##
select.rows <- function(x, k) {
  y <- x[seq(nrow(ppt)-(k-1),nrow(ppt)),]
  return(y)
}

source("prepare-data.R")

## standardize response variable and covariates before JAGS analysis
##
gi <- standardize(gi)
ppt <- standardize(ppt)
tmn <- standardize(tmn)

## calculate means for last decade for index calculation
##
ppt_mean <- apply(select.rows(ppt, 10), 2, mean)
tmn_mean <- apply(select.rows(tmn, 10), 2, mean)

## d.f. for student t in robust regression
##
nu <- 2

stan.data <- list(gi=gi,
                  year=year,
                  indiv=indiv,
                  ppt=ppt,
                  tmn=tmn,
                  ppt_mean=ppt_mean,
                  tmn_mean=tmn_mean,
                  n_obs=n.obs,
                  n_years=n.years,
                  n_indiv=n.indiv,
                  n_sites=n.sites,
                  n_months=n.months,
                  nu=nu)
stan.pars <- c("beta_0",
               "beta_ppt",
               "beta_tmn",
               "mu_year",
               "mu_indiv",
               "mu_site",
               "var_resid",
               "var_indiv",
               "var_site",
               "idx_site",
               "log_lik")
fit <- stan(file="growth-increment-with-site.stan",
            data=stan.data,
            pars=stan.pars,
            iter=n.iter,
            warmup=n.burnin,
            thin=n.thin,
            chains=n.chains,
            cores=n.cores,
            control=list(adapt_delta=0.95,
              max_treedepth=20))
opt.old <- options(width=120)
sink("results-stan.txt", append=TRUE)
cat("\n\n\n\n",
    "With site effect (and inverse gamma priors...)\n",
    "**********************************************\n",
    sep="")
print(fit,
      pars=c("beta_ppt", "beta_tmn", "var_resid", "var_indiv", "var_site"),
      digits_summary=3)
sink()
options(opt.old)

