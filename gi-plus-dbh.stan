data {
  // shared
  //
  int<lower=0> n_sites;

  // for growth increment component
  //
  int<lower=0> n_years;
  int<lower=0> n_indiv;
  int<lower=0> n_months;
  matrix[n_years,n_months] ppt;
  matrix[n_years,n_months] tmn;
  matrix[n_indiv,n_years] gi;
  int<lower=0> site_gi[n_indiv];

  // for dbh component
  //
  int<lower=0> n_obs;
  int<lower=0> n_species;
  vector[n_obs] dbh_2;
  vector[n_obs] dbh_1;
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
  vector[n_months] beta_ppt;
  vector[n_months] beta_tmn;
  vector[n_indiv] mu_indiv;
  vector[n_sites] mu_site;
  real<lower=0> sigma_indiv;
  real<lower=0> sigma_site_gi;
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;

  // for dbh component
  //
  real beta_0_dbh;
  real beta_size;
  real beta_height_ratio;
  real gamma_radiation;
  real gamma_slope;
  real gamma_aspect;
  real gamma_twi;
  real<lower=0> sigma_resid;
  real<lower=0> sigma_site_dbh;
  real<lower=0> sigma_species;
  vector[n_sites] eps_site;
  vector[n_species] eps_species;
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
  vector[n_obs] log_mu_dbh_inc;
  // prior mean for intercept for growth increment
  // derived from dbh component
  //
  vector[n_sites] beta_0_gi;

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
  for (i in 1:n_obs) {
    log_mu_dbh_inc[i] <- beta_size*tree_size[i]
                         + beta_height_ratio*height_ratio[i]
                         + gamma_radiation*radiation[site_dbh[i]]
                         + gamma_slope*slope[site_dbh[i]]
                         + gamma_aspect*aspect[site_dbh[i]]
                         + gamma_twi*twi[site_dbh[i]]
                         + eps_species[species[i]]
                         + eps_site[site_dbh[i]];
  }
  mu_indiv_dbh <- dbh_1 + exp(log_mu_dbh_inc);

  // prior mean of intercept for growth increment model derived
  // from site effect of dbh model
  //
  // divided by n_years to scale as annual increment
  //
  beta_0_gi <- exp(eps_site)/n_years;
}
model {
  // priors for growth increment component
  //
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
  beta_size ~ normal(0.0, 1.0);
  beta_height_ratio ~ normal(0.0, 1.0);
  gamma_radiation ~ normal(0.0, 1.0);
  gamma_slope ~ normal(0.0, 1.0);
  gamma_aspect ~ normal(0.0, 1.0);
  gamma_twi ~ normal(0.0, 1.0);
  sigma_resid ~ normal(0.0, 1.0);
  sigma_site_dbh ~ normal(0.0, 1.0);
  sigma_species ~ normal(0.0, 1.0);

  // likelihood for growth increment component
  //
  for (i in 1:n_indiv) {
    mu_indiv[i] ~ normal(mu_site[site_gi[i]], sigma_indiv);
  }
  for (i in 1:n_sites) {
    mu_site[i] ~ normal(beta_0_gi[i], sigma_site_gi);
  }
  // individual site x year combinations
  //
  for (i in 1:n_indiv) {
    gi[i] ~ multi_normal(mu_year_indiv[i], Sigma);
  }

  // likelihood for dbh component
  //
  eps_species ~ normal(0.0, sigma_species);
  eps_site ~ normal(beta_0_dbh, sigma_site_dbh);
  dbh_2 ~ normal(mu_indiv_dbh, sigma_resid);
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
    log_lik_dbh[i] <- normal_log(mu_indiv_dbh[i], sigma_resid)
  }
}
