library(rstan)

rm(list=ls())

stan.data <- list()
stan.pars <- c("omega")
fit <- stan(file="check-lkj.stan",
            data=stan.data,
            pars=stan.pars,
            iter=2500,
            warmup=1250,
            thin=1,
            chains=4,
            cores=4)
print(fit, digits=3)

