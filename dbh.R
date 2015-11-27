library(R2jags)
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

## set up data vectors for JAGS analysis
##
dbh.1 <- standardize(dbh$T1_BasalArea)
dbh.2 <- standardize(dbh$T2_BasalArea)
dbh.inc <- standardize(dbh$T2_BasalArea - dbh$T1_BasalArea)
size <- standardize(dbh$Tree.height)
height.ratio <- standardize(dbh$height.ratio)
plot <- as.numeric(dbh$plot)
species <- as.numeric(dbh$Species)
radiation <- standardize(plot.level$radiation)
slope <- standardize(plot.level$slope)
aspect <- standardize(plot.level$aspect)
twi <- standardize(plot.level$twi)
n.obs <- nrow(dbh)
n.plots <- length(unique(dbh$plot))
stopifnot(n.plots == max(plot))
n.species <- length(unique(dbh$Species))
stopifnot(n.species == max(species))

if (measurement.error) {
  jags.data <- c("dbh.1",
                 "dbh.2",
                 "size",
                 "height.ratio",
                 "plot",
                 "species",
                 "radiation",
                 "slope",
                 "aspect",
                 "twi",
                 "n.obs",
                 "n.plots",
                 "n.species")
  jags.pars <- c("beta.0",
                 "beta.size",
                 "beta.height.ratio",
                 "gamma.radiation",
                 "gamma.slope",
                 "gamma.aspect",
                 "gamma.twi",
                 "eps.plot",
                 "eps.species",
                 "var.resid",
                 "var.plot",
                 "var.species")
  model.file <- "dbh-measurement-error.jags"
} else {
  jags.data <- c("dbh.inc",
                 "size",
                 "height.ratio",
                 "plot",
                 "n.obs",
                 "n.plots")
  jags.pars <- c("beta.0",
                 "beta.size",
                 "beta.height.ratio",
                 "mu.indiv",
                 "var.resid",
                 "var.plot")
  model.file <- "dbh.jags"
}
fit <- jags(data=jags.data,
            inits=NULL,
            parameters=jags.pars,
            model.file=model.file,
            n.chains=n.chains,
            n.burnin=n.burnin,
            n.iter=n.iter,
            n.thin=n.thin,
            working.directory=".")
opt.old <- options(width=120)
print(fit, digits.summary=3)
options(opt.old)


