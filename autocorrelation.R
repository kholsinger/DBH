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

cat("Uncoupled model\n")
load("results-gi-plus-dbh-uncoupled.Rsave")
coupled <- autocorr(fit, max(year))
cat("Coupled model\n")
load("results-gi-plus-dbh.Rsave")
uncoupled <- autocorr(fit, max(year))

coupled.rho <- data.frame(Autocorrelation=coupled$rho,
                      Lag=seq(0, length(coupled$rho)-1, by=1),
                      Data="Coupled data")
uncoupled.rho <- data.frame(Autocorrelation=uncoupled$rho,
                      Lag=seq(0, length(uncoupled$rho)-1, by=1),
                      Data="Uncoupled data")
for.plot <- rbind(coupled.rho, uncoupled.rho)
for.plot <- subset(for.plot, Lag > 0)

p <- ggplot(for.plot, aes(x=Lag, y=Autocorrelation, color=Data)) +
     geom_point()
print(p)


