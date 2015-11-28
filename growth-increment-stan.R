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

## set counts
##
n.obs <- length(gi)
n.years <- end.series - start.series
n.indiv <- length(unique(indiv))

## set up table relating individual indices to sites
##
site.table <- unique(data.frame(site=gi.data$site, id=gi.data$id))
site <- as.numeric(site.table$site)
n.sites <- length(unique(site))

## prior on regression coefficients
##
tau.beta <- 1.0

## standardize response variable and covariates before JAGS analysis
##
gi <- standardize(gi)
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
                  ppt=ppt,
                  tmn=tmn,
                  ppt_mean=ppt_mean,
                  tmn_mean=tmn_mean,
                  n_obs=n.obs,
                  n_years=n.years,
                  n_indiv=n.indiv,
                  n_sites=n.sites,
                  n_months=n.months)
stan.pars <- c("beta_0",
               "beta_ppt",
               "beta_tmn",
               "mu_year",
               "mu_indiv",
               "mu_site",
               "sigma_resid",
               "sigma_indiv",
               "sigma_site",
               "idx_site")

fit <- stan(file="growth-increment.stan",
            data=stan.data,
            pars=stan.pars,
            iter=n.iter,
            warmup=n.burnin,
            thin=n.thin,
            chains=n.chains,
            cores=n.cores)
opt.old <- options(width=120)
print(fit, digits_summary=3)
options(opt.old)

save(fit, n.months, gi, year, indiv, site, start.series, end.series,
     file="results-stan.Rsave")

