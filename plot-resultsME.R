library(ggplot2)
library(reshape)

rm(list=ls())

check.residuals <- TRUE
JAGS <- FALSE

results.file <- "results-gi-plus-dbh-uncoupled.Rsave"
load(file=results.file)

## coefficient plots
##
if (JAGS) {
  library(R2jags)
  beta.0 <- fit$BUGSoutput$sims.list$beta.0
  beta.ppt <- fit$BUGSoutput$sims.list$beta.ppt
  beta.tmn <- fit$BUGSoutput$sims.list$beta.tmn
} else {
  library(rstan)
  pars <- extract(fit, pars=c("beta_0_gi", "beta_ppt", "beta_tmn"))
  beta.0 <- pars$beta_0_gi
  beta.ppt <- pars$beta_ppt
  beta.tmn <- pars$beta_tmn
}

beta.ppt <- beta.ppt[,c(seq(from=12, to=1))] # reordering of columns
beta.tmn <- beta.tmn[,c(seq(from=12, to=1))] # reordering of columns

n.reps <- nrow(beta.0)

months <- c("Jan.1", "Feb.1", "Mar.1", "Apr.1", "May.1", "Jun.1",
            "Jul.1", "Aug.1", "Sep", "Oct", "Nov", "Dec",
            "Jan", "Feb", "Mar", "Apr", "May", "Jun",
            "Jul", "Aug")
#old.months <- months[seq(from=length(months), to=length(months)-(n.months-1))]
months <- months[seq(from=length(months)-(n.months-1), to=length(months))]

par <- character(0)
month <- character(0)
values <- numeric(0)
for (i in 1:length(months)) {
  par <- c(par, rep("precipitation", n.reps))
  month <- c(month, rep(months[i], n.reps))
  values <- c(values, beta.ppt[,i])

  par <- c(par, rep("mean temperature", n.reps))
  month <- c(month, rep(months[i], n.reps))
  values <- c(values, beta.tmn[,i])
}

for.plot <- data.frame(par=par, month=month, values=values)
for.plot$month <- factor(for.plot$month, levels=months)

p <- ggplot(for.plot, aes(x=month, y=values)) +
     geom_hline(yintercept=0, linetype="dashed") +
     geom_violin(fill="blue", alpha=0.2) +
     geom_boxplot(width=0.1, fill="dark blue", outlier.colour=NA) +
     facet_grid(par ~ .)
print(p)

### other fixed effects
if (JAGS) {
  library(R2jags)
#  beta.0 <- fit$BUGSoutput$sims.list$beta.0
#  beta.ppt <- fit$BUGSoutput$sims.list$beta.ppt
#  beta.tmn <- fit$BUGSoutput$sims.list$beta.tmn
} else {
  library(rstan)
  gi.pars <- extract(fit, pars=c("beta_size_gi", "beta_size_gi_squared", 
                                 "beta_basal_area_gi", "beta_height_ratio_gi",
                                 "gamma_radiation_gi", "gamma_slope_gi",
                                 "gamma_aspect_gi", "gamma_twi_gi"))
  dbh.pars <- extract(fit, pars=c("beta_size", "beta_size_squared", 
                                  "beta_basal_area", "beta_height_ratio",
                                  "gamma_radiation_dbh", "gamma_slope_dbh",
                                  "gamma_aspect_dbh", "gamma_twi_dbh"))  
}

n.reps <- length(gi.pars$beta_size_gi)
fixed.effects <- c("size", "size^2", "plotBA", "HR", "rad", "slope",
            "asp", "wet")

submodel <- character(0)
par <- character(0)
values <- numeric(0)

for (i in 1:length(fixed.effects)) {
  submodel <- c(submodel, rep("gi", n.reps))
  par <- c(par, rep(fixed.effects[i], n.reps))
  values <- c(values, gi.pars[[i]])
  
  submodel <- c(submodel, rep("dbh", n.reps))
  par <- c(par, rep(fixed.effects[i], n.reps))
  values <- c(values, dbh.pars[[i]])
}

for.plot <- data.frame(submodel=submodel, par=par, values=values)
for.plot$par <- factor(for.plot$par, levels=fixed.effects)

p <- ggplot(for.plot, aes(x=par, y=values)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_violin(fill="blue", alpha=0.2) +
  geom_boxplot(width=0.1, fill="dark blue", outlier.colour=NA) +
  facet_grid(submodel ~ .)
print(p)


## variation across years
##
if (JAGS) {
  mu.year <- fit$BUGSoutput$sims.list$mu.year
} else {
  mu.year <- extract(fit, "mu_year")$mu_year
}
n.reps <- nrow(mu.year)
n.yrs <- ncol(mu.year)
years <- seq(start.series+1, end.series, by=1)
year <- numeric(0)
value <- numeric(0)
for (i in 1:n.yrs) {
  year <- c(year, rep(years[i], n.reps))
  value <- c(value, mu.year[,i])
}
for.plot <- data.frame(Year=as.factor(year), Predicted=value)

p <- ggplot(for.plot, aes(x=Year, y=Predicted)) +
     geom_violin(fill="blue", alpha=0.2) +
     geom_boxplot(width=0.1, fill="dark blue", outlier.colour=NA)
print(p)

## variation across individuals
##
if (JAGS) {
  mu.indiv <- fit$BUGSoutput$sims.list$mu.indiv
} else {
  mu.indiv <- extract(fit, "mu_indiv")$mu_indiv
}
n.reps <- nrow(mu.indiv)
n.inds <- ncol(mu.indiv)
indivs <- seq(1, n.inds, by=1)
site.id <- numeric(0)
indiv <- numeric(0)
value <- numeric(0)
for (i in 1:n.inds) {
  site.id <- c(site.id, rep(site[i], n.reps))
  indiv <- c(indiv, rep(i, n.reps))
  value <- c(value, mu.indiv[,i])
}
for.plot <- data.frame(Site=as.factor(site.id),
                       Individual=as.factor(indiv),
                       Predicted=value)

p <- ggplot(for.plot, aes(x=Individual, y=Predicted)) +
     geom_violin(fill="blue", alpha=0.2) +
     geom_boxplot(width=0.1, fill="dark blue", outlier.colour=NA) +
     facet_wrap(~ Site, scales="free_x")
print(p)

## variation across sites
##
if (JAGS) {
  mu.site <- fit$BUGSoutput$sims.list$mu.site
} else {
  mu.site <- extract(fit, "mu_site")$mu_site
}
n.reps <- nrow(mu.site)
n.sites <- ncol(mu.site)
sites <- seq(1, n.inds, by=1)
site <- numeric(0)
value <- numeric(0)
for (i in 2:n.sites) {
  site <- c(site, rep(i, n.reps))
  value <- c(value, mu.site[,i])
}
for.plot <- data.frame(Site=as.factor(site), Predicted=value)

p <- ggplot(for.plot, aes(x=Site, y=Predicted)) +
     geom_violin(fill="blue", alpha=0.2) +
     geom_boxplot(width=0.1, fill="dark blue", outlier.colour=NA)
print(p)


## residual plots (using posterior mean)
##
## reload data since year gets used above
##
if (check.residuals) {
  load(file=results.file)
  if (JAGS) {
    mu.year.indiv <- fit$BUGSoutput$mean$mu.year.indiv
  } else {
    mu.year.indiv.sims <- extract(fit, c("mu_year_indiv"))$mu_year_indiv
    dims <- dim(mu.year.indiv.sims)
    mu.year.indiv <- matrix(nrow=dims[2], ncol=dims[3])
    for (i in 1:dims[2]) {
      for (j in 1:dims[3]) {
        mu.year.indiv[i,j] <- mean(mu.year.indiv.sims[,i,j])
      }
    }
    gi.cens <- apply(extract(fit, pars=c("gi_cens"))$gi_cens, 2, mean)
  }
  n.obs <- length(gi)
  residual.ts <- matrix(nrow=dims[2], ncol=dims[3])
  ## uncensored observations
  ##
  Observed <- gi
  Predicted <- numeric(n.obs)
  Residual <- numeric(n.obs)
  for (i in 1:n.obs) {
    Predicted[i] <- mu.year.indiv[year[i], indiv[i]]
    residual.ts[year[i],indiv[i]] <- Observed[i] - Predicted[i]
  }
  ## censored observations
  ##
  n.obs.cens <- length(gi.cens)
  Observed <- c(Observed, rep(lower_bound, n.obs.cens))
  Predicted <- c(Predicted, numeric(n.obs.cens))
  Residual <- c(Residual, numeric(n.obs.cens))
  for (i in 1:n.obs.cens) {
    Predicted[i+n.obs] <- mu.year.indiv[year.cens[i], indiv.cens[i]]
    residual.ts[year.cens[i],indiv.cens[i]] <- Observed[i+n.obs] - Predicted[i+n.obs]
  }
  ## Residual
  ##
  Residual <- Observed - Predicted

  for.plot <- data.frame(Observed=Observed, Predicted=Predicted,
                         Residual=Residual)
  for.plot <- subset(for.plot, !is.na(Observed))

  p <- ggplot(for.plot, aes(x=Observed, y=Predicted)) +
       geom_abline(intercept=0, slope=1, linetype="dashed") +
       geom_point()
  print(p)

  p <- ggplot(for.plot, aes(x=Predicted, y=Residual)) +
       geom_hline(yintercept=0, linetype="dashed") +
       geom_point()
  print(p)
  ## Residual time series
  ##
  indivs <- seq(1, n.inds, by=1)
  site.id <- numeric(0)
  indiv <- numeric(0)
  year <- numeric(0)
  value <- numeric(0)
  n.yrs <- dims[2]
  yrs <- seq(start.series+1, end.series, by=1)
  lag.max <- 15
  lag <- numeric(0)
  ind <- numeric(0)
  auto.corr <- numeric(0)
  for (i in 1:n.inds) {
    acf.tmp <- acf(residual.ts[,i], na.action=na.pass, lag.max=lag.max)$acf[,,1]
    lag <- c(lag, seq(1, length(acf.tmp), by=1))
    ind <- c(ind, rep(i, length(acf.tmp), by=1))
    auto.corr <- c(auto.corr, acf.tmp)
  }
  for.plot <- data.frame(Lag=lag, Individual=ind, ACF=auto.corr)
  p <- ggplot(for.plot, aes(x=ACF)) + geom_histogram(binwidth=0.05) +
       geom_vline(xintercept=0, lty=2) + facet_wrap(~ Lag)
  print(p)
}
