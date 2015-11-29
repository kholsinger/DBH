library(reshape)
library(R2jags)

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
  set.seed(1)
}

## select last k rows from x
##
select.rows <- function(x, k) {
  y <- x[seq(nrow(ppt)-(k-1),nrow(ppt)),]
  return(y)
}

source("prepare-data.R")

## standardize response variable and covariates before JAGS analysis
##
gi <- standardize(gi)
ppt <- standardize(ppt)
tmn <- standardize(tmn)

## calculate means for last decade for index calculation
##
ppt.mean <- apply(select.rows(ppt, 10), 2, mean)
tmn.mean <- apply(select.rows(tmn, 10), 2, mean)

## prior for precision gon regression coefficients
##
tau.beta <- 1.0

## d.f. for student t in robust regression
##
nu <- 2

jags.data <- c("gi",
               "year",
               "indiv",
               "site",
               "ppt",
               "tmn",
               "ppt.mean",
               "tmn.mean",
               "n.obs",
               "n.years",
               "n.indiv",
               "n.sites",
               "n.months",
               "tau.beta",
               "nu")
jags.pars <- c("beta.0",
               "beta.ppt",
               "beta.tmn",
               "mu.year",
               "mu.indiv",
               "mu.site",
               "var.resid",
               "var.indiv",
               "var.site",
               "idx.site")
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
     file="results.Rsave")
