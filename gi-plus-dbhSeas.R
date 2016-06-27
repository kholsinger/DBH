library(reshape)
library(rstan)

rm(list=ls())

debug <- FALSE
compare <- FALSE
uncoupled <- FALSE
coupled <- TRUE
write.results.file <- TRUE
save <- TRUE

save <- save & !debug
write.results.file <- write.results.file & !debug

base_year <- 2004 # note that base_year is also defined on line 280

if (coupled) {
  model.file <- "gi-plus-dbhM2.stan"
} else if (uncoupled) {
  model.file <- "gi-plus-dbh-uncoupledM2.stan"
} else {
  stop("Mis-specification of model")
}

## MCMC settings
##
n.burnin <- 1250
n.iter <- 2500
n.thin <- 1
n.chains <- 4
if (debug) {
  n.burnin <- 50
  n.iter <- 100
  n.chains <- 2
}

rstan_options(auto_write = TRUE)
n.cores <- min(n.chains, parallel::detectCores())

## select last k rows from x
##
select.rows <- function(x, k) {
  y <- x[seq(nrow(ppt)-(k-1),nrow(ppt)),]
  return(y)
}


is.na.in.row <- function(x) {
  sum(is.na(x)) > 0
}
is.number.in.row <- function(x) {
  sum(!is.na(x)) > 0
}
not.is.na.in.row <- function(x) {
  !is.na.in.row(x)
}
strip.rwl <- function(x) {
  y <- sub(".rwl", "", x)
  return(y)
}
extract.tree.number <- function(x, site) {
  y <- sub(site, "",x)
  y <- sub("[abABN]", "", y, perl=TRUE)
  y <- as.numeric(substring(y, 1, 4))
  return(y)
}

## calculate R^2 (Gelman et al. 2006)
##
calc.r2 <- function(fit, y, n.years, n.indiv, n.iter, verbose=FALSE) {
  mu.year <- extract(fit, pars=c("mu_year"))$mu_year
  mu.indiv <- extract(fit, pars=c("mu_indiv"))$mu_indiv
  v.theta <- numeric(n.iter)
  v.eps <- numeric(n.iter)
  ## R^2
  ##
  for (iter in 1:n.iter) {
    if (verbose) {
      if ((iter %% 100) == 0) {
        cat(".")
        flush(stdout())
      }
      if (iter == n.iter) {
        cat(iter, "\n")
        flush(stdout())
      }
    }
    mu.year.indiv <- matrix(nrow=n.indiv, ncol=n.years)
    mu <- numeric(n.indiv*n.years)
    eps <- numeric(n.indiv*n.years)
    for (i in 1:n.indiv) {
      for (j in 1:n.years) {
        mu.year.indiv[i,j] <- mu.indiv[iter,i] + mu.year[iter,j]
      }
    }
    for (i in 1:n.indiv) {
      for (j in 1:n.years) {
        mu[i] <- mu.year.indiv[i,j]
        eps[i] <- y[i,j] - mu[i]
      }
    }
    v.theta[iter] <- var(mu+eps)
    v.eps[iter] <- var(eps)
  }
  r2 <- 1 - mean(v.eps)/mean(v.theta)
  ## lambda
  ##
  e.v.eps <- mean(v.eps)
  for (iter in 1:n.iter) {
    if (verbose) {
      if ((iter %% 100) == 0) {
        cat(".")
        flush(stdout())
      }
      if (iter == n.iter) {
        cat(iter, "\n")
        flush(stdout())
      }
    }
    for (i in 1:n.indiv) {
      for (j in 1:n.years) {
        mu.year.indiv[i,j] <- mean(mu.indiv[,i] + mu.year[,j])
      }
    }
    for (i in 1:n.indiv) {
      for (j in 1:n.years) {
        mu[i] <- mu.year.indiv[i,j]
        eps[i] <- y[i,j] - mu[i]
      }
    }
  }
  v.e.eps <- var(eps)
  lambda <- 1 - v.e.eps/e.v.eps
  return(list(r2=r2, lambda=lambda))
}

to.basal.area <- function(dbh) {
  area <- pi*((dbh/2.0)^2)
  return(area)
}

initial.basal.area <- function(obs, gi, base) {## obs is T1_DBH, gi.raw are ring widths, base is the T1 year
  ## calculating BA backwards from year T1 (2004) to start.series+1 (1982)
  ## initial basal area thus refers to 1981 basal area
  n.trees <- nrow(gi)
  pred <- numeric(n.trees)
  for (i in 1:n.trees) {
    dbh <- obs[i] ## each tree's T1_DBH
    for (yr in base[i]:(start.series+1)) {## from each tree's T1 year (2004) to 1981+1
      ## gi is radial growth increment
      ##
      dbh <- dbh - 2.0*gi[i, as.character(yr)]*0.1
    }
    pred[i] <- to.basal.area(dbh)
  }
  return(pred)
}

check.initial.basal.area <- function(obs, gi, base) {
  n.trees <- nrow(gi)
  pred <- numeric(n.trees)
  pred.dbh <- numeric(n.trees)
  for (i in 1:n.trees) {
    dbh <- obs[i]
    for (yr in base[i]:(start.series+1)) {
      ## gi is radial growth increment
      ##
      dbh <- dbh - 2.0*gi[i, as.character(yr)]*0.1
    }
    pred[i] <- to.basal.area(dbh)
    pred.dbh[i] <- dbh
  }
  return(data.frame(initial=pred,
                    initial.dbh=pred.dbh))
}

get.size.series <- function(obs, gi.raw, base, start, end) {
  tree.size <- initial.basal.area(obs, gi.raw, base) ## basal area of each tree in 1981
  ## gi will contain growth *area* increments as calculated from
  ## radial growth increments
  ##
  gi <- matrix(nrow=nrow(gi.raw), ncol=ncol(gi.raw))
  colnames(gi) <- colnames(gi.raw)
  current.size <- matrix(nrow=nrow(gi.raw), ncol=ncol(gi.raw)) ## basal area of each tree in each year
  colnames(current.size) <- colnames(gi.raw)
  n.trees <- nrow(gi.raw)
  for (i in 1:n.trees) {
    old.size <- tree.size[i] ## basal area in 1981
    ## convert size (area) to radius
    ##
    current.radius <- sqrt(old.size/pi) ## this is 1981 radius
    for (yr in (start+1):end) { ## from 1982 to 2013
      current.radius <- current.radius + gi.raw[i, as.character(yr)]*0.1
      ## to.basal.area takes diameter as argument, not radius
      ##
      new.size <- to.basal.area(2.0*current.radius)
      gi[i, as.character(yr)] <- new.size - old.size ## basal area increment calculation
      current.size[i, as.character(yr)] <- new.size ## basal area for each tree, each year
      old.size <- new.size
    }
    cat(i, " ", as.character(base[i]), ": ")
    cat(current.size[i, as.character(base[i])], " ")
    cat(to.basal.area(obs[i]), "\n")
    stopifnot(abs(current.size[i, as.character(base[i])]
                  - to.basal.area(obs[i])) < 1.0e-10)
  }
  return(list(gi=gi,
              current.size=current.size))
}

basal.area.inc <- function(dbh.1, dbh.2) {
  initial <- to.basal.area(dbh.1)
  final <- to.basal.area(dbh.2)
  inc <- final - initial
  return(inc)
}

source("prepare-data.R")
source("dbh-process-data.R")

## exclude plot 154 for the time being
##
dbh <- subset(dbh, plot!="154")
data <- subset(data, site!="154.rwl")

## exclude single outlier in dbh data
dbh <- subset(dbh, Tree.number!="691")

## exclude single outlier in gi data
data <- subset(data, id!="1560562a")

## standardize now so that these covariates match between the gi data
## and the dbh data
##
dbh$Tree.height <- standardize(dbh$Tree.height)
dbh$T1_BasalArea <- standardize(dbh$T1_BasalArea)
dbh$height.ratio <- standardize(dbh$height.ratio)
dbh$total <- standardize(dbh$radiation)
dbh$Slope <- standardize(dbh$Slope)
dbh$Aspect <- standardize(dbh$Aspect)
dbh$TWI <- standardize(dbh$SagaTWI)
dbh$BA2004 <- standardize(dbh$BA2004)
## exclude any individuals with NA for row
##
data <- data[apply(data, 1, not.is.na.in.row),]
## merge into single data frame for growth increment analysis,
## but first have to add plot and Tree.number columns to match
##
data$plot <- strip.rwl(data$site)
data$Tree.number <- 0
for (i in 1:nrow(data)) {
  data$Tree.number[i] <- extract.tree.number(data$id[i], data$plot[i])
}
before <- nrow(data)
data <- merge(data, dbh)
after <- nrow(data)
## make sure levels of each factor reflect the dropped plot
##
dbh <- droplevels(dbh)
data <- droplevels(data)
## drop outliers
##
dbh <- subset(dbh, !(plot==158 & Tree.number==669))
data <- subset(data, !(plot==137 & Tree.number==2175))

##
## observations
##

## growth increment data
##
## exclude 2014 from gi data
##
base_year <- dbh$T1year
gi.years <- years[1:(length(years)-1)]
gi.raw <- data[,gi.years]
site_gi <- as.numeric(as.factor(data$plot))
size.series <- get.size.series(data$T1_DBH, gi.raw, base_year,
                               start.series, end.series)
## This funny line is here because we were experimenting with analyses using
## increments in basal area (size.series$gi) instead of increments in radius.
## Variation in basal area increments is primarily among trees rather than within
## trees, and the among tree patterns don't make sense.
##
gi <- gi.raw
basal_area_gi <- data$BA2004
tree_size_gi <- data$T1_BasalArea
height_ratio_gi <- data$height.ratio
radiation_gi <- data$total
slope_gi <- data$Slope
aspect_gi <- data$Aspect
twi_gi <- data$TWI
n_sites_gi <- length(unique(site_gi))

n.indiv <- nrow(data)

### seasonal climate variables
ppt.data <- data.frame(ppt)
tmn.data <- data.frame(tmn)
tmn.warm <- tmn.data$X1 + tmn.data$X2 + tmn.data$X9 + tmn.data$X10 + tmn.data$X11
ppt.JFM <- ppt.data$X5 + ppt.data$X6 + ppt.data$X7
ppt.cool <- ppt.data$X3 + ppt.data$X4 + ppt.data$X5 + ppt.data$X6 + ppt.data$X7
ppt.warm <- ppt.data$X1 + ppt.data$X2 + ppt.data$X9 + ppt.data$X10 + ppt.data$X11

ppt <- standardize(ppt)
tmn <- standardize(tmn)
tmn.warm <- standardize(tmn.warm)
ppt.JFM <- standardize(ppt.JFM)
ppt.cool <- standardize(ppt.cool)
ppt.warm <- standardize(ppt.warm)



## dbh data
##
dbh_inc <- standardize(dbh$DBH_inc)
basal_area <- dbh$BA2004
tree_size <- dbh$T1_BasalArea
height_ratio <- dbh$height.ratio
site_dbh <- as.numeric(as.factor(dbh$plot))
species <- as.numeric(dbh$Species)
radiation <- dbh$total
slope <- dbh$Slope
aspect <- dbh$Aspect
twi <- dbh$TWI
n_sites_dbh <- length(unique(site_dbh))
n_obs <- nrow(dbh)
n_species <- length(unique(dbh$Species))
stopifnot(n_species == max(species))

stan.data <- list(gi=gi,
                  site_gi=site_gi,
                  basal_area_gi=basal_area_gi,
                  tree_size_gi=tree_size_gi,
                  height_ratio_gi=height_ratio_gi,
                  radiation_gi=radiation_gi,
                  slope_gi=slope_gi,
                  aspect_gi=aspect_gi,
                  twi_gi=twi_gi,
                  ppt_cool=ppt.cool,
                  tmn_warm=tmn.warm,
                  n_years=n.years,
                  n_indiv=n.indiv,
                  n_sites_gi=n_sites_gi,
                  n_months=n.months,
                  dbh_inc=dbh_inc,
                  basal_area=basal_area,
                  tree_size-tree_size,
                  height_ratio=height_ratio,
                  slope=slope,
                  aspect=aspect,
                  twi=twi,
                  site_dbh=site_dbh,
                  species=species,
                  n_sites_dbh=n_sites_dbh,
                  n_obs=n_obs,
                  n_species=n_species)
stan.pars <- c("beta_0_gi",
               "beta_ppt_cool",
               "beta_tmn",
               "beta_tmn_size",
               "mu_year",
               "mu_indiv",
               "mu_site",
               "sigma_indiv",
               "sigma_site_gi",
               "eta_sq",
               "rho_sq",
               "sigma_sq",
               "beta_0_dbh",
               "beta_size",
               "beta_size_squared",
               "beta_basal_area",
               "beta_height_ratio",
               "gamma_radiation_dbh",
               "gamma_slope_dbh",
               "gamma_aspect_dbh",
               "gamma_twi_dbh",
               "beta_size_gi",
               "beta_size_gi_squared",
               "beta_basal_area_gi",
               "beta_height_ratio_gi",
               "gamma_radiation_gi",
               "gamma_slope_gi",
               "gamma_aspect_gi",
               "gamma_twi_gi",
               "sigma_resid",
               "sigma_site_dbh",
               "sigma_species",
               "log_lik_gi",
               "log_lik_dbh")
if (coupled) {
  stan.pars <- c(stan.pars,
                 "alpha")
}
fit <- stan(file=model.file,
            data=stan.data,
            pars=stan.pars,
            iter=n.iter,
            warmup=n.burnin,
            thin=n.thin,
            chains=n.chains,
            cores=n.cores,
            adapt.delta=0.95)
opt.old <- options(width=120)
if (write.results.file) {
  sink("results-gi-plus-dbh.txt", append=TRUE, split=TRUE)
}
cat("\n\n************************************************************\n",
    sep="")
cat(date(), "\n",
    "With Gaussian process model for individuals\n",
    sep="")
if (coupled) {
  cat("   gi regression proportional to dbh regression\n", sep="")
} else {
  cat("   gi and dbh regression uncoupled\n", sep="")
}
if (use.detrended) {
  cat("     using detrended data\n", sep="")
} else {
  cat("     using raw data\n", sep="")
}
cat(before - after,
    " individuals discarded from growth increment analysis\n",
    "because of missing covariate data (plot, radiation, slope,\n",
    "aspect, twi)\n",
    sep="")
print.pars <- c("mu_site",
                "sigma_indiv",
                "sigma_site_gi",
                "eta_sq",
                "rho_sq",
                "sigma_sq",
                "beta_0_dbh",
                "beta_size",
                "beta_size_squared",
                "beta_basal_area",
                "beta_height_ratio",
                "gamma_radiation_dbh",
                "gamma_slope_dbh",
                "gamma_aspect_dbh",
                "gamma_twi_dbh",
                "beta_0_gi",
                "beta_ppt_cool",
                "beta_tmn",
                "beta_tmn_size",
                "beta_size_gi",
                "beta_size_gi_squared",
                "beta_basal_area_gi",
                "beta_height_ratio_gi",
                "gamma_radiation_gi",
                "gamma_slope_gi",
                "gamma_aspect_gi",
                "gamma_twi_gi",
                "sigma_resid",
                "sigma_site_dbh",
                "sigma_species")
if (coupled) {
  print.pars <- c(print.pars,
                  "alpha")
}
print(fit,
      pars=print.pars,
      digits_summary=3)
##
## R^2
##
## Not working for now
##
if (0) {
  r2 <- calc.r2(fit, gi, n.years, n.indiv, n.iter)
  cat("\n\n",
      "R^2:    ", round(r2$r2, 3), "\n",
      "lambda: ", round(r2$lambda, 3), sep="")
}
 if (write.results.file) {
  sink()
}
options(opt.old)
if (save) {
  if (coupled) {
    save(fit, data, ppt, tmn, dbh, n.months, start.series, end.series,
         file="results-gi-plus-dbh.Rsave")
  } else {
    save(fit, n.months, gi, year, indiv, site, start.series, end.series,
         file="results-gi-plus-dbh-uncoupled.Rsave")
  }
}

if (compare) {
  stan.data <- list(n_obs=n_obs,
                    n_plots=n_sites_dbh,
                    n_species=n_species,
                    dbh_inc=dbh_inc,
                    tree_size=tree_size,
                    height_ratio=height_ratio,
                    radiation=radiation,
                    slope=slope,
                    aspect=aspect,
                    twi=twi,
                    plot=site_dbh,
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
                "sigma_species",
                "eps_plot",
                "eps_species")
  fit <- stan(file="dbh.stan",
              data=stan.data,
              pars=stan.par,
              iter=n.iter,
              warmup=n.burnin,
              thin=n.thin,
              chains=n.chains,
              cores=n.cores)
  opt.old <- options(width=120)
  cat("\n\nResults from dbh.stan\n")
  print(fit, digits_summary=3)
  options(opt.old)

  stan.data <- list(gi=gi,
                    site=site_gi,
                    ppt=ppt,
                    tmn=tmn,
                    n_years=n.years,
                    n_indiv=n.indiv,
                    n_sites=n.sites,
                    n_months=n.months)
  stan.pars <- c("beta_0",
                 "beta_ppt",
                 "beta_tmn",
                 "mu_year",
                 "mu_indiv",
                 "sigma_indiv",
                 "sigma_site",
                 "eta_sq",
                 "rho_sq",
                 "sigma_sq",
                 "log_lik")
  fit <- stan(file="growth-increment-with-site.stan",
              data=stan.data,
              pars=stan.pars,
              iter=n.iter,
              warmup=n.burnin,
              thin=n.thin,
              chains=n.chains,
              cores=n.cores)
  opt.old <- options(width=120)
  cat("\n\nResults from growth-increment-with-site.stan\n")
  print(fit,
        pars=c("beta_0",
               "beta_ppt",
               "beta_tmn",
               "sigma_indiv",
               "sigma_site",
               "eta_sq",
               "rho_sq",
               "sigma_sq"),
        digits_summary=3)
  options(opt.old)
}

