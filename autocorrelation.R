library(rstan)
library(ggplot2)

rm(list=ls())

autocorr <- function(fit, n.years) {
  pars <- extract(fit, pars=c("eta_sq", "rho_sq", "sigma_sq"))
  eta_sq <- pars$eta_sq
  rho_sq <- pars$rho_sq
  sigma_sq <- pars$sigma_sq
  n.samp <- length(eta_sq)
  sigma <- matrix(nrow=n.samp, ncol=n.years+1)
  rho <- matrix(nrow=n.samp, ncol=n.years+1)
  for (i in 1:n.samp) {
    sigma[i,1] <- eta_sq[i] + sigma_sq[i]
    rho[i,1] <- 1.0
    for (j in 2:(n.years+1)) {
      sigma[i,j] <- eta_sq[i]*exp(-rho_sq[i]*((j-1)^2))
      rho[i,j] <- sigma[i,j]/sigma[i,1]
    }
  }
  cat("Autocovariance...\n")
  for (i in 1:(n.years+1)) {
    output <- sprintf("  lag %2d: %5.3f (%5.3f, %5.3f)\n",
                      i-1,
                      mean(sigma[,i]),
                      quantile(sigma[,i], probs=0.025),
                      quantile(sigma[,i], probs=0.975))
    cat(output)
  }
  cat("Autocorrelation...\n")
  for (i in 1:(n.years+1)) {
    output <- sprintf("  lag %2d: %5.3f (%5.3f, %5.3f)\n",
                      i-1,
                      mean(rho[,i]),
                      quantile(rho[,i], probs=0.025),
                      quantile(rho[,i], probs=0.975))
    cat(output)
  }
  return(list(sigma=apply(sigma, 2, mean), rho=apply(rho, 2, mean)))
}

cat("With raw data\n")
load("results-stan-raw-data.Rsave")
raw <- autocorr(fit, max(year))
cat("With detrended data\n")
load("results-stan.Rsave")
detrend <- autocorr(fit, max(year))

raw.rho <- data.frame(Autocorrelation=raw$rho,
                      Lag=seq(0, length(raw$rho)-1, by=1),
                      Data="Raw data")
det.rho <- data.frame(Autocorrelation=detrend$rho,
                      Lag=seq(0, length(detrend$rho)-1, by=1),
                      Data="Detrended data")
for.plot <- rbind(raw.rho, det.rho)
for.plot <- subset(for.plot, Lag > 0)

p <- ggplot(for.plot, aes(x=Lag, y=Autocorrelation, color=Data)) +
     geom_point()
print(p)


