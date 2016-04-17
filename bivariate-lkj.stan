data {
  int<lower=0> n_cov;
  int<lower=0> n_obs;
  matrix[n_obs,n_cov] covars;
}
parameters {
  vector[n_cov] mu;
  vector<lower=0>[n_cov] tau;
  corr_matrix[n_cov] rho;
//  cholesky_factor_corr[n_cov] rhoC;
}
transformed parameters {
  matrix[n_cov,n_cov] Sigma;
//  matrix[n_cov,n_cov] rho;

//  rho <- multiply_lower_tri_self_transpose(rhoC);
  for (i in 1:n_cov) {
    for (j in i:n_cov) {
      Sigma[i,j] <- rho[i,j]*tau[i]*tau[j];
    }
  }
}
model {
  tau ~ cauchy(0.0, 2.5);
  rho ~ lkj_corr(2.0);
//  rhoC ~ lkj_corr_cholesky(2.0);
  for (i in 1:n_obs) {
    covars[i] ~ multi_normal(mu, quad_form_diag(rho, tau));
//    covars[i] ~ multi_normal_cholesky(mu, rhoC);
  }
}
