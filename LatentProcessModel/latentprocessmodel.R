#### code for running JAGS model fusing tree-ring and DBH data

### load libraries

### set MCMC settings

### source files that prepare data

### set indexing limits

### priors for regression coefficients

### standardize predictors

### define jags data

### define jags pars

### call fit

fit <- jags(data=
            inits=
            parameters=
            model.file=
            n.chains=
            n.burnin=
            n.iter=
            n.thin=
              )

save(fit, ...file="LPMresults.Rsave")

