library(rstan)

rm(list=ls())

debug <- FALSE
wishart <- FALSE

## nrows <- 100
## ncols <- 4
## covars <- matrix(nrow=nrows, ncol=ncols)
## for (i in 1:nrows) {
##   covars[i, ] <- rnorm(ncols)
## }

## print(kappa(covars))

## load covariates and convert to format used in JAGS analysis
##
source("prepare-data.R")

months <- c("Jan.1", "Feb.1", "Mar.1", "Apr.1", "May.1", "Jun.1",
            "Jul.1", "Aug.1", "Sep.1", "Oct.1", "Nov.1", "Dec.1",
            "Jan.2", "Feb.2", "Mar.2", "Apr.2", "May.2", "Jun.2",
            "Jul.2", "Aug.2")
months <- months[seq(from=length(months), to=length(months)-(n.months-1))]

covars <- cbind(ppt, tmn)
colnames(covars) <- c(paste("PPT.", months, sep=""),
                      paste("TmeanC.", months, sep=""))

covars <- standardize(covars)
covars <- covars[,1:10]

n_obs <- nrow(covars)
n_cov <- ncol(covars)

W <- diag(x=1, nrow=n_cov, ncol=n_cov)
nu <- n_cov + 6
if (wishart) {
  stan.data <- list(n_cov=n_cov,
                    n_obs=n_obs,
                    covars=covars,
                    W=W,
                    nu=nu)
  stan.pars <- c("rho",
                 "Sigma",
                 "mu")
  fit <- stan(file="bivariate.stan",
              data=stan.data,
              pars=stan.pars,
              iter=2500,
              warmup=1250,
              thin=1,
              chains=4,
              cores=4)
} else {
  stan.data <- list(n_cov=n_cov,
                    n_obs=n_obs,
                    covars=covars)
  stan.pars <- c("rho",
                 "Sigma",
                 "mu")
  fit <- stan(file="bivariate-lkj.stan",
              data=stan.data,
              pars=stan.pars,
              iter=2500,
              warmup=1250,
              thin=1,
              chains=4,
              cores=4)
}
corr <- extract(fit, pars=c("rho"))$rho
opt.old <- options(width=180)
print(fit, digits=3, pars="rho")
options(opt.old)
