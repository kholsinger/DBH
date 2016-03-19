library(ggplot2)
library(plyr)

rm(list=ls())

## 80% credible intervals
##
alpha <- 0.2

lo <- function(x) {
  p <- quantile(x, alpha/2.0)
  return(p)
}

hi <- function(x) {
  p <- quantile(x, 1.0-alpha/2.0)
  return(p)
}

summarize <- function(x) {
  tmp <- ddply(x, c("Fire", "Substrate"), summarise,
               mean=mean(value),
               lo=lo(value),
               hi=hi(value))
  print(tmp, digits=3)
}

load("results-gi-plus-dbh-multi-correlated-with-size.Rsave")
mu.site <- extract(fit, pars=c("mu_site"))$mu_site
n.sites <- ncol(mu.site)
n.samps <- nrow(mu.site)

key <- read.csv("plot-key.csv", header=TRUE)

Fire <- character(0)
Substrate <- character(0)
value <- numeric(0)
for (i in 1:n.sites) {
  if (key$fire[i] == 1) {
    label <- "Fire"
  } else if (key$fire[i] == 0) {
    label <- "No Fire"
  } else {
    stop(paste("Error in", key$fire[i], "\n", sep=" "))
  }
  Fire <- c(Fire, rep(label, n.samps))
  if (key$substrate[i] == "T") {
    label <- "Tuff"
  } else if (key$substrate[i] =="P") {
    label <- "Pumice"
  } else if (key$substrate[i] == "A") {
    label <- "Alluvium"
  } else {
    stop(paste("Error in", key$substrate[i], "\n", sep=" "))
  }
  Substrate <- c(Substrate, rep(label, n.samps))
  value <- c(value, mu.site[,i])
}

for.plot <- data.frame(Fire=Fire,
                       Substrate=Substrate,
                       value=value)
p <- ggplot(for.plot, aes(x=value)) +
  geom_density(fill="blue", alpha=0.2) +
  xlab("Plot random effect") +
  facet_grid(Fire ~ Substrate)
print(p)

p <- ggplot(for.plot, aes(x=Substrate, y=value)) +
  geom_violin(fill="blue", alpha=0.2) +
  geom_boxplot(width=0.1, fill="dark blue", outlier.colour=NA) +
  stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) +
  ylab("Plot random effect") +
  facet_grid(Fire ~ .)
print(p)

p <- ggplot(for.plot, aes(x=Substrate, y=value)) +
  geom_violin(fill="blue", alpha=0.2) +
  geom_boxplot(width=0.1, fill="dark blue", outlier.colour=NA) +
  stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) +
  ylab("Plot random effect") +
  facet_grid(. ~ Fire)
print(p)
