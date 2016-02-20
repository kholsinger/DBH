data {
  int<lower=0> n_years;
  int<lower=0> n_indiv;
  int<lower=0> n_sites;
  int<lower=0> n_months;
  matrix[n_years,n_months] ppt;
  matrix[n_years,n_months] tmn;
  matrix[n_indiv,n_years] gi;
  int<lower=0> site[n_indiv];
  // for index calculations
  vector[n_months] ppt_mean;
  vector[n_months] tmn_mean;
}
parameters {
  vector[n_months] beta_ppt;
  vector[n_months] beta_tmn;
  vector[n_indiv] mu_indiv;
  vector[n_sites] mu_site;
  real beta_0;
  real<lower=0> sigma_indiv;
  real<lower=0> sigma_site;
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;
}
transformed parameters {
  matrix[n_indiv,n_years] mu_year_indiv;
  vector[n_years] mu_year;
  real<lower=0> rho_sq;
  cov_matrix[n_years] Sigma;

  rho_sq <- inv(inv_rho_sq);
  // beta_0 incorporated into intercept for mu_year_indiv through mu_indiv
  //
  mu_year <- ppt*beta_ppt + tmn*beta_tmn;
  // expectation for individual j in year i is sum of year
  // and indivdidual effects
  //
  for (i in 1:n_indiv) {
    for (j in 1:n_years) {
      mu_year_indiv[i,j] <- mu_indiv[i] + mu_year[j];
     }
  }
  // covariance matrix for Gaussian process
  //
  for (i in 1:(n_years-1)) {
    for (j in (i+1):n_years) {
      Sigma[i,j] <- eta_sq*exp(-rho_sq*pow(i-j,2));
      Sigma[j,i] <- Sigma[i,j];
    }
  }
  for (i in 1:n_years) {
    Sigma[i,i] <- eta_sq + sigma_sq;
  }
}
model {
  // priors
  //
  beta_0 ~ normal(0.0, 1.0);
  beta_ppt ~ normal(0.0, 1.0);
  beta_tmn ~ normal(0.0, 1.0);
  sigma_indiv ~ normal(0.0, 1.0);
  sigma_site ~ normal(0.0, 1.0);
  eta_sq ~ normal(0.0, 1.0);
  inv_rho_sq ~ normal(0.0, 1.0);
  sigma_sq ~ normal(0.0, 1.0);

  // likelihood
  //
  for (i in 1:n_indiv) {
    mu_indiv[i] ~ normal(mu_site[site[i]], sigma_indiv);
  }
  mu_site ~ normal(beta_0, sigma_site);
  // individual site x year combinations
  //
  for (i in 1:n_indiv) {
    gi[i] ~ multi_normal(mu_year_indiv[i], Sigma);
  }
}
generated quantities {
  vector[n_indiv] log_lik;
//   vector[n_sites] idx_site;

  // beta_0 incorporated through prior mean on mu_site
  //
//   idx_site <- ppt_mean'*beta_ppt + tmn_mean'*beta_tmn + mu_site;
  // calculate log likelihoods
  //
  for (i in 1:n_indiv) {
    log_lik[i] <- multi_normal_log(gi[i],
                                   mu_year_indiv[i],
                                   Sigma);
  }
}
