## construct bivariate associations between each covariate and growth
## increment from partial regression coefficients and covariate correlations
##
## FIX ME: Iterations hard coded in call to JAGS. Must match iterations
## available from growth-increment.R. Need to store JAGS parameters in
## .Rsave file and retrieve them to set equivalent parameters here
##
library(ggplot2)
library(R2jags)
library(rstan)

rm(list=ls())

debug <- FALSE

clean <- function(x) {
  y <- sub(".1|.2", "", x)
  return(y)
}

## load posterior estimates from Stan analysis
results.file <- "results-gi-plus-dbh-uncoupled.Rsave"
load(file=results.file)
beta.ppt <- extract(fit, pars=c("beta_ppt"))$beta_ppt
beta.tmn <- extract(fit, pars=c("beta_tmn"))$beta_tmn
n.iter <- nrow(beta.ppt)

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
Omega <- diag(x=1.0, nrow=ncol(covars), ncol=ncol(covars))
nu <- nrow(Omega) + 2

n.obs <- nrow(covars)
n.cov <- ncol(covars)
tau <- 1.0
jags.data <- c("covars",
               "n.obs",
               "n.cov",
               "tau",
               "Omega",
               "nu")
jags.par <- c("rho")
fit <- jags(data=jags.data,
            inits=NULL,
            parameters=jags.par,
            model.file="bivariate.jags",
            n.chains=5,
            n.burnin=1000,
            n.iter=2000,
            n.thin=1,
            working.directory=".")
corr <- fit$BUGSoutput$sims.list$rho
opt.old <- options(width=180)
print(fit, digits.summary=3)
options(opt.old)

covariate <- character(0)
month <- character(0)
correlation <- numeric(0)
for (i in 1:ncol(corr)) {
  if (i <= length(months)) {
    covar.name <- "PPT"
    month.name <- clean(months[i])
   } else {
    covar.name <- "TmeanC"
    month.name <- clean(months[i-length(months)])
  }
  for (k in 1:n.iter) {
    tmp <- 0.0
    for (j in 1:ncol(corr)) {
      if (j <= length(months)) {
        part.corr <- beta.ppt[k,j]
      } else {
        part.corr <- beta.tmn[k,(j-length(months))]
      }
      tmp <- tmp + corr[k,i,j]*part.corr
    }
    correlation <- c(correlation, tmp)
  }
  covariate <- c(covariate, rep(covar.name, n.iter))
  month <- c(month, rep(month.name, n.iter))
}

for.plot <- data.frame(Month=month,
                       Covariate=covariate,
                       Correlation=correlation)
for.plot$Month <- factor(for.plot$Month, levels=clean(months))

p <- ggplot(for.plot, aes(x=Month, y=Correlation)) +
     geom_hline(yintercept=0, linetype="dashed") +
     geom_boxplot(width=0.2, fill="dark blue", outlier.colour="dark blue") +
     stat_summary(fun.y=median, geom="point", fill="white",
                  shape=21, size=2.5) +
     facet_grid(Covariate ~ .)
print(p)

