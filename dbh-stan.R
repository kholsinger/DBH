library(rstan)
library(plyr)

rm(list=ls())

debug <- FALSE
measurement.error <- TRUE

## MCMC settings
##
n.burnin <- 1000
n.iter <- 2000
n.thin <- 1
n.chains <- 5
if (measurement.error) {
  n.burnin <- 5000
  n.iter <- 30000
  n.thin <- 25
}
if (debug) {
  n.chains <- 2
  ## to allow replication of results across runs in JAGS
  ##
  set.seed(1)
}

source("dbh-process-data.R")

## set up data vectors for Stan analysis
##
dbh_1 <- standardize(dbh$T1_BasalArea)
dbh_2 <- standardize(dbh$T2_BasalArea)
dbh_inc <- standardize(dbh$T2_BasalArea - dbh$T1_BasalArea)
tree_size <- standardize(dbh$Tree.height)
height_ratio <- standardize(dbh$height.ratio)
plot <- as.numeric(dbh$plot)
species <- as.numeric(dbh$Species)
radiation <- standardize(plot.level$radiation)
slope <- standardize(plot.level$slope)
aspect <- standardize(plot.level$aspect)
twi <- standardize(plot.level$twi)
n_obs <- nrow(dbh)
n_plots <- length(unique(dbh$plot))
stopifnot(n_plots == max(plot))
n_species <- length(unique(dbh$Species))
stopifnot(n_species == max(species))

stan.data <- list(n_obs=n_obs,
                  n_plots=n_plots,
                  n_species=n_species,
                  dbh_1=dbh_1,
                  dbh_2=dbh_2,
                  tree_size=tree_size,
                  height_ratio=height_ratio,
                  radiation=radiation,
                  slope=slope,
                  aspect=aspect,
                  twi=twi,
                  plot=plot,
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
            chains=n.chains)
opt.old <- options(width=120)
print(fit, digits_summary=3)
options(opt.old)


