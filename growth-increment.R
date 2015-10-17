library(dplR)
library(bootRes)
library(pspline)
library(reshape)
library(R2jags)

rm(list=ls())

debug <- FALSE

## load detr() and helper functions
source("detr.R")
source("allfunctionsdavid.v13.R")

month.no <- function(year, month, start.year) {
  num <- (year-start.year)*12 + month
  return(num)
}

## N.B.: VPD series does not include 2014. Requires manual adjustments below
##
start.series <- 1981
end.series <- 2014
## final.month: last month of data relevant to current year
##   growth increment
final.month <- 8
## n.months: number of months included as covariates
##
n.months <- 8

## MCMC settings
##
n.burnin <- 1000
n.iter <- 6000
n.thin <- 1
if (debug) {
  n.chains <- 2
  ## to allow replication of results across runs in JAGS
  ##
  set.seed(1)
} else {
  n.chains <- 5
}

## set paths
plot.data.path <- "plot-data"

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
  data.tr.det[[i]] <- window(data.tr.det[[i]],
                             start=start.series+1,
                             end=end.series)
}

## climate data
##
## end.series re-defined to deal with VPD mismatch
##
end.series <- 2013
prism <- read.table("PRISM_MCN.txt", header = T)      #starts Jan 1981
prism <- subset(prism, year <= end.series)
vpd <- ts(read.table("vpd_MCC_cru1901_1913.txt", header = T),
          start = 1901,
          frequency = 12)
vpd81 <- window(vpd, start = start.series, end = end.series+1, extend = TRUE)
##
data.climate <- cbind(prism, head(as.vector(vpd81), -1))
colnames(data.climate) <- c(colnames(prism), "vpd")
data.climate <- subset(data.climate, !is.na(vpd))
rm(prism, vpd, vpd81)

## translate data into format required for JAGS analysis
##

## individual growth increment data
## only for first set of plot data at the moment
##
data <- data.frame(t(data.tr.det[[1]]))
colnames(data) <- seq(from=start.series+1, to=end.series+1, by=1)
data$"2014" <- 2014
gi.data <- reshape(data,
                   idvar="id",
                   ids=row.names(data),
                   times=as.numeric(colnames(data)),
                   timevar="yr",
                   varying=list(colnames(data)),
                   direction="long",
                   v.names="gi")
gi.data <- subset(gi.data, yr <= end.series)
gi.data$yr <- gi.data$yr - start.series
## exclude growth increments with estimate=0
##
gi.data <- subset(gi.data, gi > 0)
##
## extract for JAGS
##
gi <- gi.data$gi
year <- gi.data$yr
indiv <- as.numeric(as.factor(gi.data$id))

## weather data by year
##
ppt <- matrix(nrow=end.series-start.series, ncol=n.months)
tmn <- matrix(nrow=end.series-start.series, ncol=n.months)
tmx <- matrix(nrow=end.series-start.series, ncol=n.months)
vpd <- matrix(nrow=end.series-start.series, ncol=n.months)
for (i in (start.series+1):end.series) {
  yr <- i - start.series
  ##
  month.end <- month.no(i, final.month+1, start.series)
  for (j in 1:n.months) {
    month <- month.end-j
    if (debug) {
      cat(i, j, month, "\n")
    }
    ppt[yr,j] <- data.climate$ppt.mm[month]
    tmn[yr,j] <- data.climate$TmeanC[month]
    tmx[yr,j] <- data.climate$TmaxC[month]
    vpd[yr,j] <- data.climate$vpd[month]
  }
}
if (debug) {
  if (n.months < 9) {
    months <- c("Jan.1", "Feb.1", "Mar.1", "Apr.1", "May.1", "Jun.1",
                "Jul.1", "Aug.1", "Sep.1", "Oct.1", "Nov.1", "Dec.1",
                "Jan.2", "Feb.2", "Mar.2", "Apr.2", "May.2", "Jun.2",
                "Jul.2", "Aug.2")
    months <- months[seq(from=length(months), to=length(months)-(n.months-1))]

    rownames(ppt) <- seq(start.series+1, end.series, by=1)
    colnames(ppt) <- months
    cat(ppt["1988",],"\n")
    cat(data.climate$ppt.mm[data.climate$year==1988][final.month:(final.month+1-6)], "\n")
    cat(ppt["1998",],"\n")
    cat(data.climate$ppt.mm[data.climate$year==1998][final.month:(final.month+1-6)], "\n")
    cat(ppt["2008",],"\n")
    cat(data.climate$ppt.mm[data.climate$year==2008][final.month:(final.month+1-6)], "\n")
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
               "vpd",
               "n.obs",
               "n.years",
               "n.indiv",
               "n.months",
               "tau.beta")
jags.pars <- c("beta.0",
               "beta.ppt",
               "beta.tmn",
               "beta.vpd",
               "mu.year",
               "mu.indiv",
               "mu.year.indiv",
               "var.resid",
               "var.indiv")
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

save(fit, n.months, gi, year, indiv, start.series, end.series,
     file="results.Rsave")
