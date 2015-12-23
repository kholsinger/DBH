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
  n.burnin <- 50
  n.iter <- 100
  n.chains <- 2
  ## to allow replication of results across runs in JAGS
  ##
}

rstan_options(auto_write = TRUE)
n.cores <- min(n.chains, parallel::detectCores())

## select last k rows from x
##
select.rows <- function(x, k) {
  y <- x[seq(nrow(ppt)-(k-1),nrow(ppt)),]
  return(y)
}

## calculate R^2 (Gelman et al. 2006)
##
calc.r2 <- function(mu.year.site.sims, y, n.indiv, site, verbose=FALSE) {
  dims <- dim(mu.year.site.sims$mu_year_site)
  n.iter <- dims[1]
  n.site <- dims[2]
  n.years <- dims[3]
  v.theta <- numeric(n.iter)
  v.eps <- numeric(n.iter)
  ## R^2
  ##
  for (iter in 1:n.iter) {
    if (verbose) {
      if ((iter %% 100) == 0) {
        cat(".")
        flush(stdout())
      }
      if (iter == n.iter) {
        cat(iter, "\n")
        flush(stdout())
      }
    }
    mu.year.site <- matrix(nrow=n.site, ncol=n.years)
    mu <- numeric(n.site*n.years)
    eps <- numeric(n.site*n.years)
    for (i in 1:n.site) {
      for (j in 1:n.years) {
        mu.year.site[i,j] <- mu.year.site.sims$mu_year_site[iter,i,j]
      }
    }
    for (i in 1:n.indiv) {
      mu[i] <- mu.year.site[site[i],j]
      eps[i] <- y[i,j] - mu[i]
    }
    v.theta[iter] <- var(mu+eps)
    v.eps[iter] <- var(eps)
  }
  r2 <- 1 - mean(v.eps)/mean(v.theta)
  ## lambda
  ##
  e.v.eps <- mean(v.eps)
  for (iter in 1:n.iter) {
    if (verbose) {
      if ((iter %% 100) == 0) {
        cat(".")
        flush(stdout())
      }
      if (iter == n.iter) {
        cat(iter, "\n")
        flush(stdout())
      }
    }
    ## mu.year.site <- matrix(nrow=dims[2], ncol=dims[3])
    for (i in 1:n.site) {
      for (j in 1:n.years) {
        mu.year.site[i,j] <- mean(mu.year.site.sims$mu_year_site[,i,j])
      }
    }
    for (i in 1:n.indiv) {
      mu[i] <- mu.year.site[site[i],j]
      eps[i] <- y[i,j] - mu[i]
    }
  }
  v.e.eps <- var(eps)
  lambda <- 1 - v.e.eps/e.v.eps
  return(list(r2=r2, lambda=lambda))
}

source("prepare-data.R")

## standardize response variable and covariates before JAGS analysis
##
## extracct individual-level data
## exclude any individuals with NA for now
##
is.na.in.row <- function(x) {
  sum(is.na(x)) > 0
}
not.is.na.in.row <- function(x) {
  !is.na.in.row(x)
}
data <- data[apply(data, 1, not.is.na.in.row),]
## observations
##
gi <- data[,1:n.years]
site <- as.numeric(data$site)
n.site <- length(unique(site))

n.indiv <- nrow(gi)

ppt <- standardize(ppt)
tmn <- standardize(tmn)

## calculate means for last decade for index calculation
##
ppt_mean <- apply(select.rows(ppt, 10), 2, mean)
tmn_mean <- apply(select.rows(tmn, 10), 2, mean)

stan.data <- list(gi=gi,
                  site=site,
                  ppt=ppt,
                  tmn=tmn,
                  ppt_mean=ppt_mean,
                  tmn_mean=tmn_mean,
                  n_years=n.years,
                  n_indiv=n.indiv,
                  n_sites=n.sites,
                  n_months=n.months)
stan.pars <- c("beta_0",
               "beta_ppt",
               "beta_tmn",
               "mu_year",
               "mu_site",
               "sigma_site",
               "eta_sq",
               "rho_sq",
               "sigma_sq")
fit <- stan(file="growth-increment-with-site.stan",
            data=stan.data,
            pars=stan.pars,
            iter=n.iter,
            warmup=n.burnin,
            thin=n.thin,
            chains=n.chains,
            cores=n.cores)#,
            ## control=list(max_treedepth=20,
            ##              adapt_delta=0.95))
opt.old <- options(width=120)
sink("results-stan.txt", append=TRUE, split=TRUE)
cat("\n\n\n\n",
    "With site effect (using raw data instead of detrended)\n",
    "N.B. without threshold model for growth increments of 0\n",
    "     with Gaussian process model for individuals\n",
    "*******************************************************\n",
    sep="")
print(fit,
      pars=c("beta_0",
             "beta_ppt",
             "beta_tmn",
             "sigma_site",
             "eta_sq",
             "rho_sq",
             "sigma_sq"),
      digits_summary=3)
##
## R^2
##
## FIX ME: Have to finish fixing R^2 calculation
##
## r2 <- calc.r2(extract(fit, pars=c("mu_year_indiv")), gi, n.indiv, site)
## cat("\n\n",
##     "R^2:    ", round(r2$r2, 3), "\n",
##     "lambda: ", round(r2$lambda, 3), sep="")
sink()
options(opt.old)

if (!debug) {
  save(fit, n.months, gi, year, indiv, site, start.series, end.series,
       file="results-stan.Rsave")
}
