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
## gi <- standardize(gi)
##
## re-extracct individual-level data, recognizing censoring
##
gi.data$gi <- standardize(gi.data$gi)
gi.data$id <- as.numeric(as.factor(gi.data$id))
## lower_bound is lower bound for data
##
lower_bound <- min(gi.data$gi)
## censored observations
##
gi.data.cens <- subset(gi.data, gi==lower_bound)
gi.cens <- gi.data.cens$gi
year.cens <- gi.data.cens$yr
indiv.cens <- gi.data.cens$id
site.cens <- gi.data.cens$site
## uncensored observations
##
gi.data <- subset(gi.data, gi > lower_bound)
gi <- gi.data$gi
year <- gi.data$yr
indiv <- gi.data$id
site <- gi.data$site

n.obs <- length(gi)
n.obs.cens <- length(gi.cens)

ppt <- standardize(ppt)
tmn <- standardize(tmn)

## calculate means for last decade for index calculation
##
ppt_mean <- apply(select.rows(ppt, 10), 2, mean)
tmn_mean <- apply(select.rows(tmn, 10), 2, mean)

stan.data <- list(gi=gi,
                  year=year,
                  indiv=indiv,
                  site=site,
                  year_cens=year.cens,
                  indiv_cens=indiv.cens,
                  lower_bound=lower_bound,
                  ppt=ppt,
                  tmn=tmn,
                  ppt_mean=ppt_mean,
                  tmn_mean=tmn_mean,
                  n_obs=n.obs,
                  n_obs_cens=n.obs.cens,
                  n_years=n.years,
                  n_indiv=n.indiv,
                  n_indiv_cens=n.inidv.cens,
                  n_sites=n.sites,
                  n_months=n.months,
                  nu=nu)
stan.pars <- c("beta_0",
               "beta_ppt",
               "beta_tmn",
               "mu_year",
               "mu_indiv",
               "mu_site",
               "mu_year_indiv",
               "sigma_resid",
               "sigma_indiv",
               "sigma_site",
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
    "With site effect (using raw data instead of detrended)\n",
    "******************************************************\n",
    sep="")
print(fit,
      pars=c("beta_ppt",
             "beta_tmn",
             "sigma_resid",
             "sigma_indiv",
             "sigma_site"),
      digits_summary=3)
sink()
options(opt.old)

save(fit, n.months, gi, year, indiv, site, start.series, end.series,
     file="results-stan.Rsave")
