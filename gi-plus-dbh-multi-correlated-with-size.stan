data {
  // for growth increment component
  //
  int<lower=0> n_sites_gi;
  int<lower=0> n_years;
  int<lower=0> n_indiv;
  int<lower=0> n_months;
  matrix[n_years,n_months] ppt;
  matrix[n_years,n_months] tmn;
  matrix[n_indiv,n_years] gi;
  int<lower=0> site_gi[n_indiv];
  vector[n_indiv] tree_size_gi;
  vector[n_indiv] height_ratio_gi;
  vector[n_indiv] radiation_gi;
  vector[n_indiv] slope_gi;
  vector[n_indiv] aspect_gi;
  vector[n_indiv] twi_gi;

  // for dbh component
  //
  int<lower=0> n_sites_dbh;
  int<lower=0> n_obs;
  int<lower=0> n_species;
  vector[n_obs] dbh_inc;
  vector[n_obs] tree_size;
  vector[n_obs] height_ratio;
  vector[n_obs] radiation;
  vector[n_obs] slope;
  vector[n_obs] aspect;
  vector[n_obs] twi;
  int<lower=0> site_dbh[n_obs];
  int<lower=0> species[n_obs];
}
parameters {
  // for growth increment component
  //
  real beta_0_gi;
  vector[n_months] beta_ppt;
  vector[n_months] beta_tmn;
  vector[n_indiv] mu_indiv;
  vector[n_sites_gi] mu_site;
  real<lower=0> sigma_indiv;
  real<lower=0> sigma_site_gi;
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;

  // for dbh component
  //
  real beta_0_dbh;
  real<lower=0> sigma_resid;
  real<lower=0> sigma_site_dbh;
  real<lower=0> sigma_species;
  vector[n_sites_dbh] eps_site;
  // vector[n_species] eps_species;
  vector[n_species-1] eps_species_raw;

  // for the connection
  //
  corr_matrix[12] omega;
  vector[12] gamma;
}
transformed parameters {
  // for growth increment component of the model
  //
  matrix[n_indiv,n_years] mu_year_indiv;
  vector[n_years] mu_year;
  real<lower=0> rho_sq;
  cov_matrix[n_years] Sigma;
  // for dbh component of the model
  //
  vector[n_obs] mu_indiv_dbh;
  vector[n_obs] mu_dbh_inc;
  vector[n_obs] alpha_indiv;
  // for shared component of the model
  //
  real beta_size_dbh;
  real beta_height_ratio_dbh;
  real gamma_radiation_dbh;
  real gamma_slope_dbh;
  real gamma_aspect_dbh;
  real gamma_twi_dbh;
  real beta_size_gi;
  real beta_height_ratio_gi;
  real gamma_radiation_gi;
  real gamma_slope_gi;
  real gamma_aspect_gi;
  real gamma_twi_gi;
  vector<lower=0>[12] tau;
  vector[12] zero;
  vector[n_species] eps_species;


  for (i in 1:12) {
    zero[i] <- 0.0;
    tau[i] <- 1.0;
  }

  // for growth increment component of the model
  //
  rho_sq <- inv(inv_rho_sq);
  // beta_0_gi incorporated into intercept for mu_year_indiv through mu_indiv
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

  // for dbh component of the model
  //
  beta_size_dbh <- gamma[1];
  beta_height_ratio_dbh <- gamma[2];
  gamma_radiation_dbh <- gamma[3];
  gamma_slope_dbh <- gamma[4];
  gamma_aspect_dbh <- gamma[5];
  gamma_twi_dbh <- gamma[6];
  for (i in 1:(n_species-1)) {
    eps_species[i] <- eps_species_raw[i];
  }
  eps_species[n_species] <- -sum(eps_species_raw);
  for (i in 1:n_obs) {
    mu_dbh_inc[i] <- beta_size_dbh*tree_size[i]
                     + beta_height_ratio_dbh*height_ratio[i]
                     + gamma_radiation_dbh*radiation[i]
                     + gamma_slope_dbh*slope[i]
                     + gamma_aspect_dbh*aspect[i]
                     + gamma_twi_dbh*twi[i]
                     + eps_species[species[i]]
                     + eps_site[site_dbh[i]];
  }

  // shared component at individual level
  //
  beta_size_gi <- gamma[7];
  beta_height_ratio_gi <- gamma[8];
  gamma_radiation_gi <- gamma[9];
  gamma_slope_gi <- gamma[10];
  gamma_aspect_gi <- gamma[11];
  gamma_twi_gi <- gamma[12];
  for (i in 1:n_indiv) {
    alpha_indiv[i] <- beta_size_gi*tree_size_gi[i]
                      + beta_height_ratio_gi*height_ratio_gi[i]
                      + gamma_radiation_gi*radiation_gi[i]
                      + gamma_slope_gi*slope_gi[i]
                      + gamma_aspect_gi*aspect_gi[i]
                      + gamma_twi_gi*twi_gi[i];
  }
}
model {
  // priors for growth increment component
  //
  beta_0_gi ~ normal(0.0, 1.0);
  beta_ppt ~ normal(0.0, 1.0);
  beta_tmn ~ normal(0.0, 1.0);
  sigma_indiv ~ normal(0.0, 1.0);
  sigma_site_gi ~ normal(0.0, 1.0);
  eta_sq ~ normal(0.0, 1.0);
  inv_rho_sq ~ normal(0.0, 1.0);
  sigma_sq ~ normal(0.0, 1.0);
  // priors for dbh component
  //
  beta_0_dbh ~ normal(0.0, 1.0);
  sigma_resid ~ normal(0.0, 1.0);
  sigma_site_dbh ~ normal(0.0, 1.0);
  sigma_species ~ normal(0.0, 1.0);
  // prior for shared components
  //
  omega ~ lkj_corr(2.0);
  gamma ~ multi_normal(zero, quad_form_diag(omega, tau));

  // likelihood for growth increment component
  //
  for (i in 1:n_indiv) {
    mu_indiv[i] ~ normal(mu_site[site_gi[i]] + alpha_indiv[i], sigma_indiv);
  }
  for (i in 1:n_sites_gi) {
    mu_site[i] ~ normal(beta_0_gi, sigma_site_gi);
  }
  // individual site x year combinations
  //
  for (i in 1:n_indiv) {
    gi[i] ~ multi_normal(mu_year_indiv[i], Sigma);
  }

  // likelihood for dbh component
  //
  eps_species_raw ~ normal(0.0, sigma_species);
  eps_site ~ normal(beta_0_dbh, sigma_site_dbh);
  dbh_inc ~ normal(mu_dbh_inc, sigma_resid);
}
generated quantities {
  vector[n_indiv] log_lik_gi;
  vector[n_obs] log_lik_dbh;

  // calculate log likelihoods
  //
  for (i in 1:n_indiv) {
    log_lik_gi[i] <- multi_normal_log(gi[i],
                                      mu_year_indiv[i],
                                      Sigma);
  }
  for (i in 1:n_obs) {
    log_lik_dbh[i] <- normal_log(dbh_inc[i], mu_dbh_inc[i], sigma_resid);
  }
}
