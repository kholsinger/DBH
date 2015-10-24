library(R2jags)
library(ggplot2)
library(reshape)

rm(list=ls())

results.file <- "results.Rsave"
load(file=results.file)

## coefficient plots
##
beta.0 <- fit$BUGSoutput$sims.list$beta.0
beta.ppt <- fit$BUGSoutput$sims.list$beta.ppt
beta.tmn <- fit$BUGSoutput$sims.list$beta.tmn

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
mu.year <- fit$BUGSoutput$sims.list$mu.year
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
mu.indiv <- fit$BUGSoutput$sims.list$mu.indiv
n.reps <- nrow(mu.indiv)
n.inds <- ncol(mu.indiv)
indivs <- seq(1, n.inds, by=1)
indiv <- numeric(0)
value <- numeric(0)
for (i in 2:n.inds) {
  indiv <- c(indiv, rep(i, n.reps))
  value <- c(value, mu.indiv[,i])
}
for.plot <- data.frame(Individual=as.factor(indiv), Predicted=value)

p <- ggplot(for.plot, aes(x=Individual, y=Predicted)) +
     geom_boxplot()
print(p)

## residual plots (using posterior mean)
##
## reload data since year gets used above
##
load(file=results.file)
mu.year.indiv <- fit$BUGSoutput$mean$mu.year.indiv
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

