library(ggplot2)
library(reshape)

rm(list=ls())

check.residuals <- TRUE
JAGS <- FALSE

if (JAGS) {
  results.file <- "results.Rsave"
} else {
  results.file <- "results-stan.Rsave"
}
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
  pars <- extract(fit, pars=c("beta_0", "beta_ppt", "beta_tmn"))
  beta.0 <- pars$beta_0
  beta.ppt <- pars$beta_ppt
  beta.tmn <- pars$beta_tmn
}

n.reps <- nrow(beta.0)

months <- c("Jan.1", "Feb.1", "Mar.1", "Apr.1", "May.1", "Jun.1",
            "Jul.1", "Aug.1", "Sep.1", "Oct.1", "Nov.1", "Dec.1",
            "Jan.2", "Feb.2", "Mar.2", "Apr.2", "May.2", "Jun.2",
            "Jul.2", "Aug.2")
months <- months[seq(from=length(months), to=length(months)-(n.months-1))]

par <- character(0)
month <- character(0)
values <- numeric(0)
for (i in 1:length(months)) {
  par <- c(par, rep("PPT", n.reps))
  month <- c(month, rep(months[i], n.reps))
  values <- c(values, beta.ppt[,i])

  par <- c(par, rep("Tmean", n.reps))
  month <- c(month, rep(months[i], n.reps))
  values <- c(values, beta.tmn[,i])
}

for.plot <- data.frame(par=par, month=month, values=values)
for.plot$month <- factor(for.plot$month, levels=months)

p <- ggplot(for.plot, aes(x=month, y=values)) +
     geom_hline(yintercept=0, linetype="dashed") +
     geom_boxplot() +
     facet_grid(par ~ .)
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
     geom_boxplot()
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
     geom_boxplot() +
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
     geom_boxplot()
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
  }
  n.obs <- length(gi)
  Observed <- gi
  Predicted <- numeric(n.obs)
  Residual <- numeric(n.obs)
  for (i in 1:n.obs) {
    Predicted[i] <- mu.year.indiv[year[i], indiv[i]]
  }
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
}
