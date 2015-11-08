library(R2jags)
library(plyr)

rm(list=ls())

debug <- FALSE

## MCMC settings
##
n.burnin <- 1000
n.iter <- 2000
n.thin <- 1
if (debug) {
  n.chains <- 2
  ## to allow replication of results across runs in JAGS
  ##
  set.seed(1)
} else {
  n.chains <- 5
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
dbh <- read.csv("mCNallDBH.csv", header=TRUE)
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
##
## get height ratio by plot and species
## warnings turned off because of empty plot/species combinations
##
ow <- options("warn" = -1)
dbh.sum <- ddply(dbh, c("plot","Species"),
                 summarise,
                 max=max(Tree.height, na.rm=TRUE))
options(ow)
dbh <- merge(dbh, dbh.sum)
dbh$height.ratio <- dbh$Tree.height/dbh$max
rm(dbh.sum)
##
## exclude rows lacking T2_DBH or Tree.height
##
dbh <- subset(dbh, !is.na(T1_DBH) & !is.na(T2_DBH) & !is.na(Tree.height))
dbh <- droplevels(dbh)

## set up data vectors for JAGS analysis
##
dbh.inc <- dbh$T2_DBH - dbh$T1_DBH
size <- standardize(dbh$T1_DBH)
height.ratio <- standardize(dbh$height.ratio)
plot <- as.numeric(dbh$plot)
species <- as.numeric(dbh$Species)
n.obs <- nrow(dbh)
n.plots <- length(unique(dbh$plot))
stopifnot(n.plots == max(plot))
n.species <- length(unique(dbh$Species))
stopifnot(n.species == max(species))

jags.data <- c("dbh.inc",
               "size",
               "height.ratio",
               "plot",
               "species",
               "n.obs",
               "n.plots",
               "n.species")
jags.pars <- c("beta.0",
               "beta.size",
               "beta.height.ratio",
               "var.resid",
               "var.plot",
               "var.species")
fit <- jags(data=jags.data,
            inits=NULL,
            parameters=jags.pars,
            model.file="dbh.jags",
            n.chains=n.chains,
            n.burnin=n.burnin,
            n.iter=n.iter,
            n.thin=n.thin,
            working.directory=".")
opt.old <- options(width=120)
print(fit, digits.summary=3)
options(opt.old)


