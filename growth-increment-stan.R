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

## calculate R^2 (Gelman et al. 2006)
##
calc.r2 <- function(mu.year.indiv.sims, y, y.cens, verbose=FALSE) {
  dims <- dim(mu.year.indiv.sims$mu_year_indiv)
  v.theta <- numeric(dims[1])
  v.eps <- numeric(dims[1])
  ## R^2
  ##
  for (iter in 1:dims[1]) {
    if (verbose) {
      if ((iter %% 100) == 0) {
        cat(".")
        flush(stdout())
      }
      if (iter == dims[1]) {
        cat(iter, "\n")
        flush(stdout())
      }
    }
    mu.year.indiv <- matrix(nrow=dims[2], ncol=dims[3])
    for (i in 1:dims[2]) {
      for (j in 1:dims[3]) {
        mu.year.indiv[i,j] <- mu.year.indiv.sims$mu_year_indiv[iter,i,j]
      }
    }
    mu <- numeric(n.obs+n.obs.cens)
    eps <- numeric(n.obs+n.obs.cens)
    for (i in 1:n.obs) {
      mu[i] <- mu.year.indiv[year[i], indiv[i]]
      eps[i] <- y[i] - mu[i]
    }
    for (i in 1:n.obs.cens) {
      mu[i+n.obs] <- mu.year.indiv[year.cens[i], indiv.cens[i]]
      eps[i+n.obs] <- y.cens[i] - mu[i+n.obs]
    }
    v.theta[iter] <- var(mu+eps)
    v.eps[iter] <- var(eps)
  }
  r2 <- 1 - mean(v.eps)/mean(v.theta)
  ## lambda
  ##
  e.v.eps <- mean(v.eps)
  for (iter in 1:dims[1]) {
    if (verbose) {
      if ((iter %% 100) == 0) {
        cat(".")
        flush(stdout())
      }
      if (iter == dims[1]) {
        cat(iter, "\n")
        flush(stdout())
      }
    }
    mu.year.indiv <- matrix(nrow=dims[2], ncol=dims[3])
    for (i in 1:dims[2]) {
      for (j in 1:dims[3]) {
        mu.year.indiv[i,j] <- mean(mu.year.indiv.sims$mu_year_indiv[,i,j])
      }
    }
    mu <- numeric(n.obs+n.obs.cens)
    eps <- numeric(n.obs+n.obs.cens)
    for (i in 1:n.obs) {
      mu[i] <- mu.year.indiv[year[i], indiv[i]]
      eps[i] <- y[i] - mu[i]
    }
  }
  v.e.eps <- var(eps)
  lambda <- 1 - v.e.eps/e.v.eps
  return(list(r2=r2, lambda=lambda))
}

source("prepare-data.R")

## standardize response variable and covariates before JAGS analysis
##
## re-extracct individual-level data, recognizing censoring
##
gi.data$gi <- standardize(sqrt(gi.data$gi))
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
## uncensored observations
##
gi.data <- subset(gi.data, gi > lower_bound)
gi <- gi.data$gi
year <- gi.data$yr
indiv <- gi.data$id

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
                  n_sites=n.sites,
                  n_months=n.months)
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
               "log_lik",
               "gi_cens")
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
sink("results-stan.txt", append=TRUE, split=TRUE)
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
##
## R^2
##
r2 <- calc.r2(extract(fit, pars=c("mu_year_indiv")), gi, rep(0, n.obs.cens))
cat("\n\n",
    "R^2:    ", round(r2$r2, 3), "\n",
    "lambda: ", round(r2$lambda, 3), sep="")
sink()
options(opt.old)

save(fit, n.months, gi, lower_bound,
     year, indiv, year.cens, indiv.cens,
     site, start.series, end.series,
     file="results-stan.Rsave")
