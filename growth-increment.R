library(reshape)
library(R2jags)

rm(list=ls())

debug <- FALSE
check.residuals <- FALSE

## MCMC settings
##
n.burnin <- 1000
n.iter <- 2000
n.thin <- 1
if (debug) {
  n.chains <- 2
  ## to allow replication of results across runs in JAGS
  ##
  set.seed(1)
} else {
  n.chains <- 5
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
site <- site.table$site
n.sites <- length(unique(site))

## prior on regression coefficients
##
tau.beta <- 1.0

## standardize response variable and covariates before JAGS analysis
##
gi <- standardize(gi)
ppt <- standardize(ppt)
tmn <- standardize(tmn)

jags.data <- c("gi",
               "year",
               "indiv",
               "site",
               "ppt",
               "tmn",
               "n.obs",
               "n.years",
               "n.indiv",
               "n.sites",
               "n.months",
               "tau.beta")
if (check.residuals) {
  jags.pars <- c("beta.0",
                 "beta.ppt",
                 "beta.tmn",
                 "mu.year",
                 "mu.indiv",
                 "mu.site",
                 "mu.year.indiv",
                 "var.resid",
                 "var.indiv",
                 "var.site")
} else {
  jags.pars <- c("beta.0",
                 "beta.ppt",
                 "beta.tmn",
                 "mu.year",
                 "mu.indiv",
                 "mu.site",
                 "var.resid",
                 "var.indiv",
                 "var.site")
}
fit <- jags(data=jags.data,
            inits=NULL,
            parameters=jags.pars,
            model.file="growth-increment.jags",
            n.chains=n.chains,
            n.burnin=n.burnin,
            n.iter=n.iter,
            n.thin=n.thin,
            working.directory=".")
opt.old <- options(width=120)
print(fit, digits.summary=3)
options(opt.old)

save(fit, n.months, gi, year, indiv, site, start.series, end.series,
     check.residuals, file="results.Rsave")
