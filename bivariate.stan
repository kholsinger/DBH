data {
  int<lower=0> n_cov;
  int<lower=0> n_obs;
  matrix[n_obs,n_cov] covars;
}
parameters {
  corr_matrix[n_cov] rho;
  vector<lower=0>[n_cov] tau;
  vector[n_cov] mu;
}
model {
  tau ~ cauchy(0.0, 2.5);
  rho ~ lkj_corr(2);
  mu ~ normal(0.0, 1.0);
  for (i in 1:n_obs) {
    covars[i] ~ multi_normal(mu, quad_form_diag(rho, tau));
  }
}
