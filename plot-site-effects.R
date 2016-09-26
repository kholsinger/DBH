library(ggplot2)
library(reshape)
library(rstan)

rm(list=ls())

results.file <- "results-gi-plus-dbh-uncoupled.Rsave"
load(file=results.file)

## variation across sites
##
mu.site <- extract(fit, "mu_site")$mu_site
n.reps <- nrow(mu.site)
n.sites <- ncol(mu.site)
site <- numeric(0)
value <- numeric(0)
for (i in 1:n.sites) {
  site <- c(site, rep(i, n.reps))
  value <- c(value, mu.site[,i])
}
for.plot <- data.frame(Site=as.factor(site),
                       Predicted=value)

p <- ggplot(for.plot, aes(x=Site, y=Predicted)) +
     geom_violin(fill="blue", alpha=0.2) +
     geom_boxplot(width=0.1, fill="dark blue", outlier.colour=NA)
print(p)
