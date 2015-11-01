## construct bivariate associations between each covariate and growth
## increment from partial regression coefficients and covariate correlations
library(ggplot2)

rm(list=ls())

debug <- FALSE

## load posterior estimates from JAGS analysis
##
results.file <- "results.Rsave"
load(file=results.file)
beta.ppt <- fit$BUGSoutput$sims.list$beta.ppt
beta.tmn <- fit$BUGSoutput$sims.list$beta.tmn
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
corr <- cor(covars)

covariate <- character(0)
month <- character(0)
correlation <- numeric(0)
for (i in 1:ncol(corr)) {
  if (i <= length(months)) {
    covar.name <- "PPT"
    month.name <- months[i]
  } else {
    covar.name <- "TmeanC"
    month.name <- months[i-length(months)]
  }
  for (k in 1:n.iter) {
    tmp <- 0.0
    for (j in 1:ncol(corr)) {
      if (j <= length(months)) {
        part.corr <- beta.ppt[k,j]
      } else {
        part.corr <- beta.tmn[k,(j-length(months))]
      }
      tmp <- tmp + corr[i,j]*part.corr
    }
    correlation <- c(correlation, tmp)
  }
  covariate <- c(covariate, rep(covar.name, n.iter))
  month <- c(month, rep(month.name, n.iter))
}

for.plot <- data.frame(Month=month,
                       Covariate=covariate,
                       Correlation=correlation)
for.plot$Month <- factor(for.plot$Month, levels=months)

p <- ggplot(for.plot, aes(x=Month, y=Correlation)) +
     geom_hline(yintercept=0, linetype="dashed") +
     geom_boxplot() +
     facet_grid(Covariate ~ .)
print(p)
