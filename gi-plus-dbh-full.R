library(reshape)
library(rstan)

rm(list=ls())

debug <- FALSE
compare <- FALSE
coupled <- FALSE
correlated <- FALSE
multi_correlated <- FALSE
multi_with_size <- TRUE
save <- TRUE

save <- save & !debug

if (multi_with_size) {
  model.file <- "gi-plus-dbh-multi-corr-wsize-full.stan"
} else if (multi_correlated) {
  model.file <- "gi-plus-dbh-multi-correlated-full.stan"
} else if (correlated) {
  model.file <- "gi-plus-dbh-correlated-full.stan"
} else if (coupled) {
  model.file <- "gi-plus-dbh-full.stan"
} else {
  model.file <- "gi-plus-dbh-uncoupled-full.stan"
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
  y <- as.numeric(substring(y, 1, 4))
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
dbh$Tree.height <- standardize(dbh$Tree.height)
dbh$height.ratio <- standardize(dbh$height.ratio)
dbh$T1_BasalArea <- standardize(dbh$T1_BasalArea)
dbh$total <- standardize(dbh$radiation)
dbh$Slope <- standardize(dbh$Slope)
dbh$Aspect <- standardize(dbh$Aspect)
dbh$Elev <- standardize(dbh$DEM)
dbh$TWI <- standardize(dbh$SagaTWI)
dbh$pBA <- standardize(dbh$BA2004)
dbh$trmi100 <- standardize(dbh$TRMI100)
dbh$trmi100 <- standardize(dbh$TRMI100)
dbh$trmi100 <- standardize(dbh$TRMI100)
dbh$trmi100 <- standardize(dbh$TRMI100)
dbh$trmi100 <- standardize(dbh$TRMI100)
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

data <- merge(data, dbh)

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
tBA_gi <- data$T1_BasalArea
height_gi <- data$Tree.height
height_ratio_gi <- data$height.ratio
radiation_gi <- data$radiation
slope_gi <- data$Slope
aspect_gi <- data$Aspect
elev_gi <- data$DEM
twi_gi <- data$TWI
fire_gi <- data$Fire2012
soil_gi <- data$soil
substrate_gi <- data$substrate
pBA_gi <- data$BA2004
trmi100_gi <- data$TRMI100
trmi250_gi <- data$TRMI250
trmi500_gi <- data$TRMI500
trmi1000_gi <- data$TRMI1000
trmi2000_gi <- data$TRMI2000
n_sites_gi <- length(unique(site_gi))

n.indiv <- nrow(data)

ppt <- standardize(ppt)
tmn <- standardize(tmn)

## dbh data
##
dbh_inc <- standardize(dbh$T2_BasalArea - dbh$T1_BasalArea)
tBA <- dbh$T1_BasalArea
height <- dbh$Tree.height
height_ratio <- dbh$height.ratio
site_dbh <- as.numeric(as.factor(dbh$plot))
species <- as.numeric(dbh$Species)
radiation <- dbh$radiation
slope <- dbh$Slope
aspect <- dbh$Aspect
elev <- dbh$DEM
twi <- dbh$SagaTWI
pBA <- dbh$BA2004
fire <- dbh$Fire2012
soil <- dbh$soil
substrate <- dbh$substrate
trmi100 <- dbh$TRMI100
trmi250 <- dbh$TRMI250
trmi500 <- dbh$TRMI500
trmi1000 <- dbh$TRMI1000
trmi2000 <- dbh$TRMI2000
n_sites_dbh <- length(unique(site_dbh))
n_obs <- nrow(dbh)
n_species <- length(unique(dbh$Species))
stopifnot(n_species == max(species))

stan.data <- list(gi=gi,
                  site_gi=site_gi,
                  tBA_gi=tBA_gi,
                  height_gi=height_gi,
                  height_ratio_gi=height_ratio_gi,
                  radiation_gi=radiation_gi,
                  slope_gi=slope_gi,
                  aspect_gi=aspect_gi,
                  elev_gi=elev_gi,
                  twi_gi=twi_gi,
                  fire_gi=fire_gi,
                  soil_gi=soil_gi,
                  substrate_gi=substrate_gi,
                  pBA_gi=pBA_gi,
                  trmi100_gi=trmi100_gi,
                  trmi250_gi=trmi250_gi,
                  trmi500_gi=trmi500_gi,
                  trmi1000_gi=trmi1000_gi,
                  trmi2000_gi=trmi2000_gi,
                  ppt=ppt,
                  tmn=tmn,
                  n_years=n.years,
                  n_indiv=n.indiv,
                  n_sites_gi=n_sites_gi,
                  n_months=n.months,
                  dbh_inc=dbh_inc,
                  tBA=tBA,
                  height=height,
                  height_ratio=height_ratio,
                  slope=slope,
                  aspect=aspect,
                  elev=elev,
                  twi=twi,
                  fire=fire,
                  soil=soil,
                  substrate=substrate,
                  pBA=pBA,
                  trmi100=trmi100,
                  trmi250=trmi250,
                  trmi500=trmi500,
                  trmi1000=trmi1000,
                  trmi2000=trmi2000,
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
               "beta_tba",
               "beta_ht",
               "beta_height_ratio",
               "gamma_radiation_dbh",
               "gamma_slope_dbh",
               "gamma_aspect_dbh",
               "gamma_elev_dbh",
               "gamma_twi_dbh",
               "gamma_fire_dbh",
               "gamma_soil_dbh",
               "gamma_substrate_dbh",
               "gamma_pba_dbh",
               "gamma_trmi100_dbh",
               "gamma_trmi250_dbh",
               "gamma_trmi500_dbh", 
               "gamma_trmi1000_dbh",
               "gamma_trmi2000_dbh",
               "gamma_radiation_gi",
               "gamma_slope_gi",
               "gamma_aspect_gi",
               "gamma_elev_gi",
               "gamma_twi_gi",
               "gamma_fire_gi",
               "gamma_soil_gi",
               "gamma_substrate_gi",
               "gamma_pba_gi",
               "gamma_trmi100_gi",
               "gamma_trmi250_gi",
               "gamma_trmi500_gi", 
               "gamma_trmi1000_gi",
               "gamma_trmi2000_gi",
               "sigma_resid",
               "sigma_site_dbh",
               "sigma_species",
               "log_lik_gi",
               "log_lik_dbh")
if (multi_correlated) {
  stan.pars <- c(stan.pars, "omega")
} else if (correlated) {
  stan.pars <- c(stan.pars, "omega_radiation",
                 "omega_slope",
                 "omega_aspect",
                 "omega_elev",
                 "omega_twi",
                 "omega_fire",
                 "omega_soil",
                 "omega_substrate",
                 "omega_pba",
                 "omega_trmi100",
                 "omega_trmi250",
                 "omega_trmi500",
                 "omega_trmi1000",
                 "omega_trmi2000")
} else if (coupled) {
  stan.pars <- c(stan.pars, "alpha")
}
if (multi_with_size) {
  stan.data <- list(gi=gi,
                    site_gi=site_gi,
                    tBA_gi=tBA_gi,
                    height_gi=height_gi,
                    height_ratio_gi=height_ratio_gi,
                    radiation_gi=radiation_gi,
                    slope_gi=slope_gi,
                    aspect_gi=aspect_gi,
                    elev_gi=elev_gi,
                    twi_gi=twi_gi,
                    fire_gi=fire_gi,
                    soil_gi=soil_gi,
                    substrate_gi=substrate_gi,
                    pBA_gi=pBA_gi,
                    trmi100_gi=trmi100_gi,
                    trmi250_gi=trmi250_gi,
                    trmi500_gi=trmi500_gi,
                    trmi1000_gi=trmi1000_gi,
                    trmi2000_gi=trmi2000_gi,
                    ppt=ppt,
                    tmn=tmn,
                    n_years=n.years,
                    n_indiv=n.indiv,
                    n_sites_gi=n_sites_gi,
                    n_months=n.months,
                    dbh_inc=dbh_inc,
                    tBA=tBA,
                    height=height,
                    height_ratio=height_ratio,
                    slope=slope,
                    aspect=aspect,
                    elev=elev,
                    twi=twi,
                    fire=fire,
                    soil=soil,
                    substrate=substrate,
                    pBA=pBA,
                    trmi100=trmi100,
                    trmi250=trmi250,
                    trmi500=trmi500,
                    trmi1000=trmi1000,
                    trmi2000=trmi2000,
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
                 "beta_tba_dbh",
                 "beta_ht_dbh",
                 "beta_height_ratio_dbh",
                 "gamma_radiation_dbh",
                 "gamma_slope_dbh",
                 "gamma_aspect_dbh",
                 "gamma_elev_dbh",
                 "gamma_twi_dbh",
                 "gamma_fire_dbh",
                 "gamma_soil_dbh",
                 "gamma_substrate_dbh",
                 "gamma_pba_dbh",
                 "gamma_trmi100_dbh",
                 "gamma_trmi250_dbh",
                 "gamma_trmi500_dbh", 
                 "gamma_trmi1000_dbh",
                 "gamma_trmi2000_dbh",
                 "beta_tba_gi",
                 "beta_ht_gi",
                 "beta_height_ratio_gi",
                 "gamma_radiation_gi",
                 "gamma_slope_gi",
                 "gamma_aspect_gi",
                 "gamma_elev_gi",
                 "gamma_twi_gi",
                 "gamma_fire_gi",
                 "gamma_soil_gi",
                 "gamma_substrate_gi",
                 "gamma_pba_gi",
                 "gamma_trmi100_gi",
                 "gamma_trmi250_gi",
                 "gamma_trmi500_gi", 
                 "gamma_trmi1000_gi",
                 "gamma_trmi2000_gi",
                 "sigma_resid",
                 "sigma_site_dbh",
                 "sigma_species",
                 "log_lik_gi",
                 "log_lik_dbh",
                 "omega")
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
cat("\n\n************************************************************\n",
    sep="")
cat(date(), "\n",
    "With Gaussian process model for individuals\n",
    sep="")
if (multi_with_size) {
  cat("   gi and dbh regression multivariate correlated, full model\n", sep="")
  cat("   tree_size and height_ratio included as covariates in gi\n", sep="")
} else if (multi_correlated) {
  cat("   gi and dbh regression multivariate correlated, full model\n", sep="")
} else if (correlated) {
  cat("   gi and dbh regression correlated, full model\n", sep="")
} else if (coupled) {
  cat("   gi regression proportional to dbh regression, full model\n", sep="")
} else {
  cat("   gi and dbh regression uncoupled, full model\n", sep="")
}
if (use.detrended) {
  cat("     using detrended data\n", sep="")
} else {
  cat("     using raw data\n", sep="")
}
#cat(before - after,
#    " individuals discarded from growth increment analysis\n",
#    "because of missing covariate data (plot, radiation, slope,\n",
#    "aspect, twi)\n",
#    sep="")
print.pars <- c("beta_0_gi",
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
                "beta_tba",
                "beta_ht",
                "beta_height_ratio",
                "gamma_radiation_dbh",
                "gamma_slope_dbh",
                "gamma_aspect_dbh",
                "gamma_elev_dbh",
                "gamma_twi_dbh",
                "gamma_fire_dbh",
                "gamma_soil_dbh",
                "gamma_substrate_dbh",
                "gamma_pba_dbh",
                "gamma_trmi100_dbh",
                "gamma_trmi250_dbh",
                "gamma_trmi500_dbh", 
                "gamma_trmi1000_dbh",
                "gamma_trmi2000_dbh",
                "gamma_radiation_gi",
                "gamma_slope_gi",
                "gamma_aspect_gi",
                "gamma_elev_gi",
                "gamma_twi_gi",
                "gamma_fire_gi",
                "gamma_soil_gi",
                "gamma_substrate_gi",
                "gamma_pba_gi",
                "gamma_trmi100_gi",
                "gamma_trmi250_gi",
                "gamma_trmi500_gi", 
                "gamma_trmi1000_gi",
                "gamma_trmi2000_gi",
                "sigma_resid",
                "sigma_site_dbh",
                "sigma_species",)
if (multi_correlated) {
  print.pars <- c(print.pars, "omega")
} else if (correlated) {
  print.pars <- c(print.pars, "omega_radiation",
                  "omega_slope",
                  "omega_aspect",
                  "omega_elev",
                  "omega_twi",
                  "omega_fire",
                  "omega_soil",
                  "omega_substrate",
                  "omega_pba",
                  "omega_trmi100",
                  "omega_trmi250",
                  "omega_trmi500",
                  "omega_trmi1000",
                  "omega_trmi2000")
} else if (coupled) {
  print.pars <- c(print.pars, "alpha")
}
if (multi_with_size) {
  print.pars <- c("beta_0_gi",
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
                  "beta_tba_dbh",
                  "beta_ht_dbh",
                  "beta_height_ratio_dbh",
                  "gamma_radiation_dbh",
                  "gamma_slope_dbh",
                  "gamma_aspect_dbh",
                  "gamma_elev_dbh",
                  "gamma_twi_dbh",
                  "gamma_fire_dbh",
                  "gamma_soil_dbh",
                  "gamma_substrate_dbh",
                  "gamma_pba_dbh",
                  "gamma_trmi100_dbh",
                  "gamma_trmi250_dbh",
                  "gamma_trmi500_dbh", 
                  "gamma_trmi1000_dbh",
                  "gamma_trmi2000_dbh",
                  "beta_tba_gi",
                  "beta_ht_gi",
                  "beta_height_ratio_gi",
                  "gamma_radiation_gi",
                  "gamma_slope_gi",
                  "gamma_aspect_gi",
                  "gamma_elev_gi",
                  "gamma_twi_gi",
                  "gamma_fire_gi",
                  "gamma_soil_gi",
                  "gamma_substrate_gi",
                  "gamma_pba_gi",
                  "gamma_trmi100_gi",
                  "gamma_trmi250_gi",
                  "gamma_trmi500_gi", 
                  "gamma_trmi1000_gi",
                  "gamma_trmi2000_gi",
                  "sigma_resid",
                  "sigma_site_dbh",
                  "sigma_species",
                  "omega")
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
if (save) {
  if (multi_with_size) {
    save(fit, data, ppt, tmn, dbh, n.months, start.series, end.series,
         file="results-gi-plus-dbh-multi-correlated-with-size.Rsave")
  } else if (multi_correlated) {
    save(fit, data, ppt, tmn, dbh, n.months, start.series, end.series,
         file="results-gi-plus-dbh-multi-correlated.Rsave")
  } else if (correlated) {
    save(fit, data, ppt, tmn, dbh, n.months, start.series, end.series,
         file="results-gi-plus-dbh-correlated.Rsave")
  } else if (coupled) {
    save(fit, data, ppt, tmn, dbh, n.months, start.series, end.series,
         file="results-gi-plus-dbh.Rsave")
  } else {
    save(fit, n.months, gi, year, indiv, site, start.series, end.series,
         file="results-gi-plus-dbh-uncoupled.Rsave")
  }
}

#### the following not updated yet...as of 3/16/16
if (compare) {
  stan.data <- list(n_obs=n_obs,
                    n_plots=n_sites_dbh,
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
                "beta_ht",
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

