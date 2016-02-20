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
source("dbh-process-data.R")

## exclude plot 154 for the time being
##
dbh <- subset(dbh, plot!="154")
data <- subset(data, site!="154.rwl")

## standardize response variable and covariates before analysis
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
strip.rwl <- function(x) {
  y <- sub(".rwl", "", x)
  return(y)
}
data <- data[apply(data, 1, not.is.na.in.row),]
## observations
##

## growth increment data
##
gi <- data[,1:n.years]
site_gi <- as.numeric(as.factor(strip.rwl(data$site)))
n.sites <- length(unique(site))

n.indiv <- nrow(gi)

ppt <- standardize(ppt)
tmn <- standardize(tmn)

## dbh data
##
dbh_1 <- standardize(dbh$T1_BasalArea)
dbh_2 <- standardize(dbh$T2_BasalArea)
dbh_inc <- standardize(dbh$T2_BasalArea - dbh$T1_BasalArea)
tree_size <- standardize(dbh$Tree.height)
height_ratio <- standardize(dbh$height.ratio)
site_dbh <- as.numeric(dbh$plot)
species <- as.numeric(dbh$Species)
radiation <- standardize(plot.level$radiation)
slope <- standardize(plot.level$slope)
aspect <- standardize(plot.level$aspect)
twi <- standardize(plot.level$twi)
n_obs <- nrow(dbh)
n_species <- length(unique(dbh$Species))
stopifnot(n_species == max(species))

stan.data <- list(gi=gi,
                  site=site,
                  ppt=ppt,
                  tmn=tmn,
                  n_years=n.years,
                  n_indiv=n.indiv,
                  n_sites=n.sites,
                  n_months=n.months,
                  dbh_1=dbh_1,
                  dbh_2=dbh_2,
                  tree_size-tree_size,
                  height_ratio=height_ratio,
                  slope=slope,
                  aspect=aspect,
                  twi=twi,
                  site_dbh=site_dbh,
                  species=species,
                  n_obs=n_obs,
                  n_species=n_species)
stan.pars <- c("beta_0_gi",
               "beta_ppt",
               "beta_tmn",
               "mu_year",
               "mu_indiv",
               "sigma_indiv",
               "sigma_site_gi",
               "eta_sq",
               "rho_sq",
               "sigma_sq",
               "beta_0_dbh",
               "beta_size",
               "beta_height_ratio",
               "gamma_radiation",
               "gamma_aspect",
               "gamma_twi",
               "sigma_resid",
               "sigma_site_dbh",
               "sigma_species",
               "log_lik_gi",
               "log_lik_dbh")
fit <- stan(file="gi-plus-dbh.stan",
            data=stan.data,
            pars=stan.pars,
            iter=n.iter,
            warmup=n.burnin,
            thin=n.thin,
            chains=n.chains,
            cores=n.cores)
opt.old <- options(width=120)
if (!debug) {
  sink("results-gi-plus-dbh.txt", append=TRUE, split=TRUE)
}
cat("************************************************************\n",
    sep="")
cat(date(), "\n",
    "With Gaussian process model for individuals\n", sep="")
if (use.detrended) {
  cat("     using detrended data\n", sep="")
} else {
  cat("     using raw data\n", sep="")
}
print(fit,
      pars=c("beta_0_gi",
             "beta_ppt",
             "beta_tmn",
             "sigma_indiv",
             "sigma_site_gi",
             "eta_sq",
             "rho_sq",
             "sigma_sq",
             "beta_0_dbh",
             "beta_size",
             "beta_height_ratio",
             "gamma_radiation",
             "gamma_aspect",
             "gamma_twi",
             "sigma_resid",
             "sigma_site_dbh",
             "sigma_species"),
      digits_summary=3)
##
## R^2
##
r2 <- calc.r2(fit, gi, n.years, n.indiv)
cat("\n\n",
    "R^2:    ", round(r2$r2, 3), "\n",
    "lambda: ", round(r2$lambda, 3), sep="")
if (!debug) {
  sink()
}
options(opt.old)

if (!debug) {
  save(fit, n.months, gi, year, indiv, site, start.series, end.series,
       file="results-gi-plus-dbh.Rsave")
}
