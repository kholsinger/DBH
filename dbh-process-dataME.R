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
dbh <- read.csv("mCNallDBH3.csv", header=TRUE)
dbh$plot <- as.factor(dbh$plot)
dbh$id <- as.factor(dbh$TreeNumber)
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
dbh <- subset(dbh, plot!="154")
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
radiation <- read.csv("mcn_plots_rsun_hours_global_radiation_utm.csv",
                      header=TRUE)
radiation$total <- apply(radiation[,grep("total", colnames(radiation))], 1, sum)

## read in and summarize DEM data
##
dem <- read.csv("mcn_plots_topographic_variables_SAGAGIS.csv")

## read in and summarize TRMI data
##
trmi <- read.csv("mcn_sample_NAD27.csv")

## merge radiation and DEM data into dbh
##
dbh <- merge(dbh, radiation, by.x="plot", by.y="ID")
dbh <- merge(dbh, dem, by.x="plot", by.y="ID")
dbh <- merge(dbh, trmi, by.x="plot", by.y="plot")

## extract plot-level covariates and index
plot.level <- data.frame(plot=dbh$plot,
                         radiation=dbh$total,
                         slope=dbh$Slope,
                         aspect=dbh$Aspect,
                         twi=dbh$TWI,
                         trmi100=dbh$TRMI100,
                         trmi250=dbh$TRMI250,
                         trmi500=dbh$TRMI500,
                         trmi1000=dbh$TRMI1000)

