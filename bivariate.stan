data {
  int<lower=0> n_cov;
  int<lower=0> n_obs;
  matrix[n_obs,n_cov] covars;
  matrix[n_cov,n_cov] W;
  real nu;
}
parameters {
  cov_matrix[n_cov] Tau;
  vector[n_cov] mu;
}
transformed parameters {
  matrix[n_cov,n_cov] rho;
  matrix[n_cov,n_cov] Sigma;

  Sigma <- inverse(Tau);
  for (i in 1:n_cov) {
    rho[i,i] <- 1.0;
    for (j in (i+1):n_cov) {
      rho[i,j] <- Sigma[i,j]/(sqrt(Sigma[i,i]*Sigma[j,j]));
      rho[j,i] <- rho[i,j];
    }
  }
}
model {
  Tau ~ wishart(nu, W);
  mu ~ normal(0.0, 1.0);
  for (i in 1:n_obs) {
    covars[i] ~ multi_normal(mu, Tau);
  }
}
