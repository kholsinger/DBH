library(R2jags)
library(plyr)

rm(list=ls())

debug <- FALSE
measurement.error <- TRUE

## MCMC settings
##
n.burnin <- 1000
n.iter <- 2000
n.thin <- 1
n.chains <- 5
if (measurement.error) {
  n.burnin <- 5000
  n.iter <- 30000
  n.thin <- 25
}
if (debug) {
  n.chains <- 2
  ## to allow replication of results across runs in JAGS
  ##
  set.seed(1)
}

standardize.vector <- function(x) {
  x.bar <- mean(x, na.rm=TRUE)
  s.d. <- sd(x, na.rm=TRUE)
  return((x-x.bar)/s.d.)
}

standardize <- function(x) {
  if (is.matrix(x)) {
    y <- matrix(nrow=nrow(x),ncol=ncol(x))
    for (i in 1:ncol(x)) {
      y[,i] <- standardize.vector(x[,i])
    }
  } else {
    y <- standardize.vector(x)
  }
  return(y)
}

## read in the DBH data and convert plot to factor
##
dbh <- read.csv("mCNallDBH2.csv", header=TRUE)
dbh$plot <- as.factor(dbh$plot)
## exclude individuals without species label
##
dbh <- subset(dbh, Species!="")
## make sure all species names are uppercase
##
dbh$Species <- factor(toupper(dbh$Species))
## exclude unknowns
##
dbh <- subset(dbh, Species!="UNK")
## correct typo in species name
##
dbh$Species[dbh$Species=="PMSE"] <- "PSME"
## analyze only PIPO, since sample sizes for others are so small
##
dbh <- subset(dbh, Species=="PIPO")
##
## get height ratio by plot and species
##
dbh.sum <- ddply(dbh, c("plot"),
                 summarise,
                 max=max(Tree.height, na.rm=TRUE))
dbh <- merge(dbh, dbh.sum)
dbh$height.ratio <- dbh$Tree.height/dbh$max
rm(dbh.sum)
##
## exclude rows lacking T2_DBH or Tree.height
##
dbh <- subset(dbh, !is.na(T1_DBH) & !is.na(T2_DBH) & !is.na(Tree.height))
dbh <- droplevels(dbh)

## read in and summarize radiation data
##
radiation <- read.csv("mcn_plots_rsun_hours_global_radiation_utm.csv",
                      header=TRUE)
radiation$total <- apply(radiation[,grep("total", colnames(radiation))], 1, sum)

## read in and summarize DEM data
##
dem <- read.csv("mcn_plots_topographic_variables_SAGAGIS.csv")

## merge radiation and DEM data into dbh
##
dbh <- merge(dbh, radiation, by.x="plot", by.y="ID")
dbh <- merge(dbh, dem, by.x="plot", by.y="ID")

## extract plot-level covariates and index
plot.level <- data.frame(plot=dbh$plot,
                         radiation=dbh$total,
                         slope=dbh$Slope,
                         aspect=dbh$Aspect,
                         twi=dbh$TWI)
plot.level <- unique(plot.level)

## set up data vectors for JAGS analysis
##
dbh.1 <- standardize(dbh$T1_DBH)
dbh.2 <- standardize(dbh$T2_DBH)
dbh.inc <- standardize(dbh$T2_DBH - dbh$T1_DBH)
size <- standardize(dbh$Tree.height)
height.ratio <- standardize(dbh$height.ratio)
plot <- as.numeric(dbh$plot)
radiation <- standardize(plot.level$radiation)
slope <- standardize(plot.level$slope)
aspect <- standardize(plot.level$aspect)
twi <- standardize(plot.level$twi)
plot.idx <- as.numeric(plot.level$plot)
n.obs <- nrow(dbh)
n.plots <- length(unique(dbh$plot))
stopifnot(n.plots == max(plot))

if (measurement.error) {
  jags.data <- c("dbh.1",
                 "dbh.2",
                 "size",
                 "height.ratio",
                 "plot",
                 "radiation",
                 "slope",
                 "aspect",
                 "twi",
                 "n.obs",
                 "n.plots")
  jags.pars <- c("beta.0",
                 "beta.size",
                 "beta.height.ratio",
                 "gamma.radiation",
                 "gamma.slope",
                 "gamma.aspect",
                 "gamma.twi",
#                 "mu.indiv",
                 "var.resid",
                 "var.plot")
  model.file <- "dbh-measurement-error.jags"
} else {
  jags.data <- c("dbh.inc",
                 "size",
                 "height.ratio",
                 "plot",
                 "n.obs",
                 "n.plots")
  jags.pars <- c("beta.0",
                 "beta.size",
                 "beta.height.ratio",
                 "mu.indiv",
                 "var.resid",
                 "var.plot")
  model.file <- "dbh.jags"
}
fit <- jags(data=jags.data,
            inits=NULL,
            parameters=jags.pars,
            model.file=model.file,
            n.chains=n.chains,
            n.burnin=n.burnin,
            n.iter=n.iter,
            n.thin=n.thin,
            working.directory=".")
opt.old <- options(width=120)
print(fit, digits.summary=3)
options(opt.old)


