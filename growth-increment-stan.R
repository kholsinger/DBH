library(reshape)
library(rstan)

rm(list=ls())

debug <- TRUE

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
calc.r2 <- function(fit, y, n.years, n.indiv, verbose=FALSE) {
  mu.year <- extract(fit, pars=c("mu_year"))$mu_year
  mu.indiv <- extract(fit, pars=c("mu_indiv"))$mu_indiv
  n.iter <- dim(mu.indiv)[1]
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
    mu.year.indiv <- matrix(nrow=n.indiv, ncol=n.years)
    mu <- numeric(n.indiv*n.years)
    eps <- numeric(n.indiv*n.years)
    for (i in 1:n.indiv) {
      for (j in 1:n.years) {
        mu.year.indiv[i,j] <- mu.indiv[iter,i] + mu.year[iter,j]
      }
    }
    for (i in 1:n.indiv) {
      for (j in 1:n.years) {
        mu[i] <- mu.year.indiv[i,j]
        eps[i] <- y[i,j] - mu[i]
      }
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
    for (i in 1:n.indiv) {
      for (j in 1:n.years) {
        mu.year.indiv[i,j] <- mean(mu.indiv[,i] + mu.year[,j])
      }
    }
    for (i in 1:n.indiv) {
      for (j in 1:n.years) {
        mu[i] <- mu.year.indiv[i,j]
        eps[i] <- y[i,j] - mu[i]
      }
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
               "mu_indiv",
               "sigma_indiv",
               "sigma_site",
               "eta_sq",
               "rho_sq",
               "sigma_sq",
               "log_lik")
fit <- stan(file="growth-increment-with-site.stan",
            data=stan.data,
            pars=stan.pars,
            iter=n.iter,
            warmup=n.burnin,
            thin=n.thin,
            chains=n.chains,
            cores=n.cores)
opt.old <- options(width=120)
sink("results-stan.txt", append=TRUE, split=TRUE)
cat("\n\n\n\n",
    "With individual effect (using raw data instead of detrended)\n",
    "N.B. without threshold model for growth increments of 0\n",
    "     with Gaussian process model for individuals\n",
    "************************************************************\n",
    sep="")
print(fit,
      pars=c("beta_0",
             "beta_ppt",
             "beta_tmn",
             "sigma_indiv",
             "sigma_site",
             "eta_sq",
             "rho_sq",
             "sigma_sq"),
      digits_summary=3)
##
## R^2
##
r2 <- calc.r2(fit, gi, n.years, n.indiv)
cat("\n\n",
    "R^2:    ", round(r2$r2, 3), "\n",
    "lambda: ", round(r2$lambda, 3), sep="")
sink()
options(opt.old)

if (!debug) {
  save(fit, n.months, gi, year, indiv, site, start.series, end.series,
       file="results-stan.Rsave")
}
