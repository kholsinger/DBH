library(dplR)
library(bootRes)
library(pspline)
library(reshape)
library(R2jags)

rm(list=ls())

## load detr() and helper functions
source("detr.R")
source("allfunctionsdavid.v13.R")

## set paths
plot.data.path <- "plot-data"

start.series <- 1981
end.series <- 2014
## time period for monthly weather data
## final.month: last month of data relevant to current year
##   growth increment
## n.months: number of months included as covariates
final.month <- 8
n.months <- 20

debug <- FALSE

if (debug) {
  n.chains <- 1
  n.burnin <- 500
  n.iter <- 1000
  n.thin <- 1
  ## to allow replication of results across runs in JAGS
  ##
  set.seed(1)
} else {
  n.chains <- 5
  n.burnin <- 5000
  n.iter <- 30000
  n.thin <- 25
}

### load data
##
## tree-ring data + detrend
##
files <- list.files(path=plot.data.path, pattern = "rwl")
data.tr <- vector(mode = "list", length = length(files))
data.tr.det <- vector(mode = "list", length = length(files))
##
for(i in 1:length(files)) {
  filename <- paste(plot.data.path, "/", files[i], sep="")
  a <- read.rwl(filename, header = FALSE)
  data.tr[[i]] <- ts(a[,1:ncol(a)],
                     start = as.numeric(rownames(a)[1]),
                     frequency = 1)
  data.tr.det[[i]] <- detr(data.tr[[i]],10)
  data.tr.det[[i]] <- window(data.tr.det[[i]], start=start.series+1, end=end.series)
}

## climate data
##
prism <- read.table("PRISM_MCN.txt", header = T)      #starts Jan 1981
vpd <- ts(read.table("vpd_MCC_cru1901_1913.txt", header = T),
          start = 1901,
          frequency = 12)
vpd81 <- window(vpd, start = start.series, end = end.series, extend = T)
##
data.climate <- cbind(prism, c(as.vector(vpd81), rep(NA,11)))
colnames(data.climate) <- c(colnames(prism), "vpd")
rm(prism, vpd, vpd81)

## translate data into format required for JAGS analysis
##

## individual growth increment data
## only for first set of plot data at the moment
##
data <- data.frame(t(data.tr.det[[1]]))
colnames(data) <- seq(from=start.series+1, to=end.series, by=1)
gi.data <- reshape(data,
                   idvar="id",
                   ids=row.names(data),
                   times=as.numeric(colnames(data)),
                   timevar="yr",
                   varying=list(colnames(data)),
                   direction="long",
                   v.names="gi")
gi.data$yr <- gi.data$yr - start.series
##
## extract for JAGS
##
gi <- gi.data$gi
year <- gi.data$yr
indiv <- as.numeric(as.factor(gi.data$id))

## weather data by year
##
ppt <- matrix(nrow=end.series-start.series, ncol=n.months+1)
tmn <- matrix(nrow=end.series-start.series, ncol=n.months+1)
tmx <- matrix(nrow=end.series-start.series, ncol=n.months+1)
vpd <- matrix(nrow=end.series-start.series, ncol=n.months+1)
for (i in (start.series+1):end.series) {
  yr <- i - start.series
  tmp <- subset(data.climate, year==(i-1) | (year==i & month <= final.month))
  ppt[yr,1] <- yr
  tmn[yr,1] <- yr
  tmx[yr,1] <- yr
  vpd[yr,1] <- yr
  for (j in 1:n.months) {
    ppt[yr,j+1] <- tmp$ppt.mm[j]
    tmn[yr,j+1] <- tmp$TmeanC[j]
    tmx[yr,j+1] <- tmp$TmaxC[j]
    vpd[yr,j+1] <- tmp$vpd[j]
  }
}

## set counts
##
n.obs <- length(gi)
n.years <- end.series - start.series
n.indiv <- length(unique(indiv))

## prior on regression coefficients
##
tau.beta <- 1.0

jags.data <- c("gi",
               "year",
               "indiv",
               "ppt",
               "tmn",
               "tmx",
               "vpd",
               "n.obs",
               "n.years",
               "n.indiv",
               "n.months",
               "tau.beta")
jags.pars <- c("beta.0",
               "beta.ppt",
               "beta.tmn",
               "beta.tmx",
               "beta.vpd",
               "tau.resid",
               "tau.indiv")
fit <- jags(data=jags.data,
            inits=NULL,
            parameters=jags.pars,
            model.file="growth-increment.jags",
            n.chains=n.chains,
            n.burnin=n.burnin,
            n.iter=n.iter,
            n.thin=n.thin,
            working.directory=".")
opt.old <- options(width=120)
print(fit, digits.summary=3)
options(opt.old)



