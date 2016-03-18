library(plyr)

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
dbh <- read.csv("MCNallDBH3.csv", header=TRUE)
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
## exclude POTR
##
dbh <- subset(dbh, Species!="POTR")
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
#radiation <- read.csv("mcn_plots_rsun_hours_global_radiation_utm.csv",
#                      header=TRUE)
#radiation$total <- apply(radiation[,grep("total", colnames(radiation))], 1, sum)

## read in and summarize DEM data
##
#dem <- read.csv("mcn_plots_topographic_variables_SAGAGIS.csv")

## read in and summarize all plot-level covariates
plot.cov <- read.csv("MCNplotCovariatesNEW.csv", header=TRUE)
plot.cov$Fire2012 <- as.factor(plot.cov$Fire2012)
plot.cov$soil <- as.factor(plot.cov$soil)
plot.cov$substrate <- as.factor(plot.cov$substrate)

## merge radiation and DEM data into dbh
##
#dbh <- merge(dbh, radiation, by.x="plot", by.y="ID")
#dbh <- merge(dbh, dem, by.x="plot", by.y="ID")
dbh <- merge(dbh, plot.cov, by.x="plot", by.y="PLOT")

## extract plot-level covariates and index
plot.level <- data.frame(plot=dbh$plot,
                         radiation=dbh$radiation,
                         slope=dbh$Slope,
                         aspect=dbh$Aspect,
                         twi=dbh$SagaTWI,
                         elev=dbh$DEM,
                         fire=dbh$Fire2012,
                         soil=dbh$soil,
                         substrate=dbh$substrate,
                         pba=dbh$BA2004,
                         TRMI100=dbh$TRMI100,
                         TRMI250=dbh$TRMI250,
                         TRMI500=dbh$TRMI500,
                         TRMI1000=dbh$TRMI1000,
                         TRMI2000=dbh$TRMI2000)
