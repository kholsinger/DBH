library(rstan)
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
## exclude plot 154 for the time being
##
## dbh <- subset(dbh, plot!="154")
## ## analyze only PIPO, since sample sizes for others are so small
## ##
## dbh <- subset(dbh, Species=="PIPO")
##
## warnings turned off because of empty plot/species combinations
##
ow <- options("warn" = -1)
## get height ratio by plot and species
##
dbh.sum <- ddply(dbh, c("plot"),
                 summarise,
                 max=max(Tree.height, na.rm=TRUE))
dbh <- merge(dbh, dbh.sum)
options(ow)
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
## plot.level <- unique(plot.level)

## set up data vectors for Stan analysis
##
dbh_1 <- standardize(dbh$T1_BasalArea)
dbh_2 <- standardize(dbh$T2_BasalArea)
dbh_inc <- standardize(dbh$T2_BasalArea - dbh$T1_BasalArea)
tree_size <- standardize(dbh$Tree.height)
height_ratio <- standardize(dbh$height.ratio)
plot <- as.numeric(dbh$plot)
species <- as.numeric(dbh$Species)
radiation <- standardize(plot.level$radiation)
slope <- standardize(plot.level$slope)
aspect <- standardize(plot.level$aspect)
twi <- standardize(plot.level$twi)
n_obs <- nrow(dbh)
n_plots <- length(unique(dbh$plot))
stopifnot(n_plots == max(plot))
n_species <- length(unique(dbh$Species))
stopifnot(n_species == max(species))

stan.data <- list(n_obs=n_obs,
                  n_plots=n_plots,
                  n_species=n_species,
                  dbh_1=dbh_1,
                  dbh_2=dbh_2,
                  tree_size=tree_size,
                  height_ratio=height_ratio,
                  radiation=radiation,
                  slope=slope,
                  aspect=aspect,
                  twi=twi,
                  plot=plot,
                  species=species)
stan.par <- c("beta_0",
              "beta_size",
              "beta_height_ratio",
              "gamma_radiation",
              "gamma_slope",
              "gamma_aspect",
              "gamma_twi",
              "sigma_resid",
              "sigma_plot",
              "sigma_species")

fit <- stan(file="dbh.stan",
            data=stan.data,
            pars=stan.par,
            chains=n.chains)
opt.old <- options(width=120)
print(fit, digits_summary=3)
options(opt.old)


