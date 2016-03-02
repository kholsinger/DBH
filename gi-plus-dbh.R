library(reshape)
library(rstan)

rm(list=ls())

debug <- FALSE
compare <- FALSE
coupled <- FALSE
correlated <- FALSE

if (correlated) {
  model.file <- "gi-plus-dbh-correlated.stan"
} else if (coupled) {
  model.file <- "gi-plus-dbh.stan"
} else {
  model.file <- "gi-plus-dbh-uncoupled.stan"
}

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


is.na.in.row <- function(x) {
  sum(is.na(x)) > 0
}
is.number.in.row <- function(x) {
  sum(!is.na(x)) > 0
}
not.is.na.in.row <- function(x) {
  !is.na.in.row(x)
}
strip.rwl <- function(x) {
  y <- sub(".rwl", "", x)
  return(y)
}
extract.tree.number <- function(x, site) {
  y <- sub(site, "",x)
  y <- sub("[abABN]", "", y, perl=TRUE)
  y <- substring(y, 1, 4)
  return(y)
}

## calculate R^2 (Gelman et al. 2006)
##
calc.r2 <- function(fit, y, n.years, n.indiv, n.iter, verbose=FALSE) {
  mu.year <- extract(fit, pars=c("mu_year"))$mu_year
  mu.indiv <- extract(fit, pars=c("mu_indiv"))$mu_indiv
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
## standardize now so that these covariates match between the gi data
## and the dbh data
##
dbh$total <- standardize(dbh$total)
dbh$Slope <- standardize(dbh$Slope)
dbh$Aspect <- standardize(dbh$Aspect)
dbh$TWI <- standardize(dbh$TWI)
## exclude any individuals with NA for row
##
data <- data[apply(data, 1, not.is.na.in.row),]
## merge into single data frame for growth increment analysis,
## but first have to add plot and Tree.number columns to match
##
data$plot <- strip.rwl(data$site)
data$Tree.number <- 0
for (i in 1:nrow(data)) {
  data$Tree.number[i] <- extract.tree.number(data$id[i], data$plot[i])
}
before <- nrow(data)
data <- merge(data, dbh)
after <- nrow(data)
## make sure levels of each factor reflect the dropped plot
##
dbh <- droplevels(dbh)
data <- droplevels(data)

##
## observations
##

## growth increment data
##
## exclude 2014 from gi data
##
gi.years <- years[1:(length(years)-1)]
gi <- data[,gi.years]
site_gi <- as.numeric(as.factor(data$plot))
radiation_gi <- data$total
slope_gi <- data$Slope
aspect_gi <- data$Aspect
twi_gi <- data$TWI
n_sites_gi <- length(unique(site_gi))

n.indiv <- nrow(data)

ppt <- standardize(ppt)
tmn <- standardize(tmn)

## dbh data
##
dbh_inc <- standardize(dbh$T2_BasalArea - dbh$T1_BasalArea)
tree_size <- standardize(dbh$Tree.height)
height_ratio <- standardize(dbh$height.ratio)
site_dbh <- as.numeric(as.factor(dbh$plot))
species <- as.numeric(dbh$Species)
radiation <- dbh$total
slope <- dbh$Slope
aspect <- dbh$Aspect
twi <- dbh$TWI
n_sites_dbh <- length(unique(site_dbh))
n_obs <- nrow(dbh)
n_species <- length(unique(dbh$Species))
stopifnot(n_species == max(species))

stan.data <- list(gi=gi,
                  site_gi=site_gi,
                  radiation_gi=radiation_gi,
                  slope_gi=slope_gi,
                  aspect_gi=aspect_gi,
                  twi_gi=twi_gi,
                  ppt=ppt,
                  tmn=tmn,
                  n_years=n.years,
                  n_indiv=n.indiv,
                  n_sites_gi=n_sites_gi,
                  n_months=n.months,
                  dbh_inc=dbh_inc,
                  tree_size-tree_size,
                  height_ratio=height_ratio,
                  slope=slope,
                  aspect=aspect,
                  twi=twi,
                  site_dbh=site_dbh,
                  species=species,
                  n_sites_dbh=n_sites_dbh,
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
               "gamma_radiation_dbh",
               "gamma_slope_dbh",
               "gamma_aspect_dbh",
               "gamma_twi_dbh",
               "gamma_radiation_gi",
               "gamma_slope_gi",
               "gamma_aspect_gi",
               "gamma_twi_gi",
               "sigma_resid",
               "sigma_site_dbh",
               "sigma_species",
               "log_lik_gi",
               "log_lik_dbh")
if (correlated) {
  stan.pars <- c(stan.pars, "omega")
} else if (coupled) {
  stan.pars <- c(stan.pars, "alpha")
}
fit <- stan(file=model.file,
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
    "With Gaussian process model for individuals\n",
    sep="")
if (correlated) {
  cat("   gi and dbh regression correlated\n", sep="")
} else if (coupled) {
  cat("   gi regression proportional to dbh regression\n", sep="")
} else {
  cat("   gi and dbh regression uncoupled\n", sep="")
}
if (use.detrended) {
  cat("     using detrended data\n", sep="")
} else {
  cat("     using raw data\n", sep="")
}
cat(before - after,
    " individuals discarded from growth increment analysis\n",
    "because of missing covariate data (plot, radiation, slope,\n",
    "aspect, twi)\n",
    sep="")
print.pars <- c("beta_0_gi",
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
                "gamma_radiation_dbh",
                "gamma_slope_dbh",
                "gamma_aspect_dbh",
                "gamma_twi_dbh",
                "gamma_radiation_gi",
                "gamma_slope_gi",
                "gamma_aspect_gi",
                "gamma_twi_gi",
                "sigma_resid",
                "sigma_site_dbh",
                "sigma_species")
if (correlated) {
  print.pars <- c(print.pars, "omega")
} else if (coupled) {
  print.pars <- c(print.pars, "alpha")
}
print(fit,
      pars=print.pars,
      digits_summary=3)
##
## R^2
##
## Not working for now
##
if (0) {
  r2 <- calc.r2(fit, gi, n.years, n.indiv, n.iter)
  cat("\n\n",
      "R^2:    ", round(r2$r2, 3), "\n",
      "lambda: ", round(r2$lambda, 3), sep="")
}
 if (!debug) {
  sink()
}
options(opt.old)
if (!debug) {
  if (correlated) {
    save(fit, data, ppt, tmn, dbh, n.months, start.series, end.series,
         file="results-gi-plus-dbh-correlated.Rsave")
  } else if (uncoupled) {
    save(fit, data, ppt, tmn, dbh, n.months, start.series, end.series,
         file="results-gi-plus-dbh-uncoupled.Rsave")
  } else {
    save(fit, n.months, gi, year, indiv, site, start.series, end.series,
         file="results-gi-plus-dbh.Rsave")
  }
}

if (compare) {
  stan.data <- list(n_obs=n_obs,
                    n_plots=n.sites,
                    n_species=n_species,
                    dbh_inc=dbh_inc,
                    tree_size=tree_size,
                    height_ratio=height_ratio,
                    radiation=radiation,
                    slope=slope,
                    aspect=aspect,
                    twi=twi,
                    plot=site_dbh,
                    species=species)
  stan.par <- c("beta_0",
                "beta_size",
                "beta_height_ratio",
                "gamma_radiation",
                "gamma_slope",
                "gamma_aspect",
                "gamma_twi",
                "sigma_resid",
                "sigma_plot",
                "sigma_species",
                "eps_plot",
                "eps_species")
  fit <- stan(file="dbh.stan",
              data=stan.data,
              pars=stan.par,
              iter=n.iter,
              warmup=n.burnin,
              thin=n.thin,
              chains=n.chains,
              cores=n.cores)
  opt.old <- options(width=120)
  cat("\n\nResults from dbh.stan\n")
  print(fit, digits_summary=3)
  options(opt.old)

  stan.data <- list(gi=gi,
                    site=site_gi,
                    ppt=ppt,
                    tmn=tmn,
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
  cat("\n\nResults from growth-increment-with-site.stan\n")
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
  options(opt.old)
}

