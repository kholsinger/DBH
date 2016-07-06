data {
  // for growth increment component
  //
  int<lower=0> n_sites_gi;
  int<lower=0> n_years;
  int<lower=0> n_indiv;
  int<lower=0> n_months;
  matrix[n_indiv,n_years] gi;
  int<lower=0> site_gi[n_indiv];
  vector[n_years] ppt_cool;
  vector[n_years] ppt_warm;
//  vector[n_years] tmn_cool;
//  vector[n_years] tmn_warm;
  vector[n_indiv] basal_area_gi;
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
  vector[n_obs] basal_area;
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
  real beta_ppt_cool;
  real beta_ppt_warm;
//  real beta_tmn_cool;
//  real beta_tmn_warm;
//  real beta_TP_warm;
//  real beta_TP_cool;
//  real beta_pcool_sq;
//  real beta_pwarm_sq;
//  real beta_tcool_sq;
//  real beta_twarm_sq;
//  real beta_pcool_htr;
//  real beta_twarm_htr;
//  real beta_tcool_htr;
//  real beta_pwarm_htr;
  real beta_pcool_BA;
//  real beta_twarm_BA;
//  real beta_tcool_BA;
  real beta_pwarm_BA;
//  real beta_twarm_size;
//  real beta_tcool_size;
  real beta_pcool_size;
  real beta_pwarm_size;
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
  real beta_size;
  real beta_size_squared;
  real beta_basal_area;
  real beta_height_ratio;
  real<lower=0> sigma_resid;
  real<lower=0> sigma_site_dbh;
  real<lower=0> sigma_species;
  vector[n_sites_dbh] eps_site;
  vector[n_species] eps_species;

  // for the connection
  //
  real<lower=0> alpha;
  real gamma_radiation_dbh;
  real gamma_slope_dbh;
  real gamma_aspect_dbh;
  real gamma_twi_dbh;
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
  // for shared component of the model, scaled by alpha from dbh component
  //
  real beta_size_gi;
  real beta_size_gi_squared;
  real beta_basal_area_gi;
  real beta_height_ratio_gi;
  real gamma_radiation_gi;
  real gamma_slope_gi;
  real gamma_aspect_gi;
  real gamma_twi_gi;

  // for growth increment component of the model
  //
  rho_sq <- inv(inv_rho_sq);
  // beta_0_gi incorporated into intercept for mu_year_indiv through mu_indiv
  //
  mu_year <- ppt_cool*beta_ppt_cool + ppt_warm*beta_ppt_warm;

//  + tmn_cool*beta_tmn_cool+ tmn_warm*beta_tmn_warm;

// interactions between T and P within a season
//+ ppt_cool*tmn_cool*beta_TP_cool + ppt_warm*tmn_warm*beta_TP_warm;

// quadratic terms climate  
//  mu_year <- ppt_cool*beta_ppt_cool + beta_pcool_sq*pow(ppt_cool, 2.0) 
//             + tmn_warm*beta_tmn_warm + beta_twarm_sq*pow(tmn_warm, 2.0) 
//             + ppt_warm*beta_ppt_warm + beta_pwarm_sq*pow(ppt_warm, 2.0) 
//             + tmn_cool*beta_tmn_cool + beta_tcool_sq*pow(tmn_cool, 2.0);             
  // expectation for individual j in year i is sum of year
  // and indivdidual effects
  //
  for (i in 1:n_indiv) {
    for (j in 1:n_years) {
      mu_year_indiv[i,j] <- mu_indiv[i] + mu_year[j] + ppt_cool[j]*basal_area_gi[i]*beta_pcool_BA + ppt_warm[j]*basal_area_gi[i]*beta_pwarm_BA + ppt_cool[j]*tree_size_gi[i]*beta_pcool_size + ppt_warm[j]*tree_size_gi[i]*beta_pwarm_size;
// climate*height ratio interactions
//+ mu_year[j] + ppt_cool[j]*height_ratio_gi[i]*beta_pcool_htr + ppt_warm[j]*height_ratio_gi[i]*beta_pwarm_htr + tmn_cool[j]*height_ratio_gi[i]*beta_tcool_htr + tmn_warm[j]*height_ratio_gi[i]*beta_twarm_htr;
// climate*plot basal area interactions
//+ ppt_cool[j]*basal_area_gi[i]*beta_pcool_BA + ppt_warm[j]*basal_area_gi[i]*beta_pwarm_BA + tmn_cool[j]*basal_area_gi[i]*beta_tcool_BA + tmn_warm[j]*basal_area_gi[i]*beta_twarm_BA;
// climate*size interactions
//        mu_year_indiv[i,j] <- mu_indiv[i] + mu_year[j] + tmn_warm[j]*tree_size_gi[i]*beta_twarm_size + tmn_cool[j]*tree_size_gi[i]*beta_tcool_size + ppt_cool[j]*tree_size_gi[i]*beta_pcool_size + ppt_warm[j]*tree_size_gi[i]*beta_pwarm_size;
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
    mu_dbh_inc[i] <- beta_size*tree_size[i]
                     + beta_size_squared*pow(tree_size[i], 2.0)
                     + beta_basal_area*basal_area[i]
                     + beta_height_ratio*height_ratio[i]
                     + gamma_radiation_dbh*radiation[i]
                     + gamma_slope_dbh*slope[i]
                     + gamma_aspect_dbh*aspect[i]
                     + gamma_twi_dbh*twi[i]
                     + eps_species[species[i]]
                     + eps_site[site_dbh[i]];
  }

  // shared component at individual level, alpha scales
  // regression coefficients from dbh component to gi component
  //
  beta_size_gi <- alpha*beta_size;
  beta_size_gi_squared <- alpha*beta_size_squared;
  beta_basal_area_gi <- alpha*beta_basal_area;
  beta_height_ratio_gi <- alpha*beta_height_ratio;
  gamma_radiation_gi <- alpha*gamma_radiation_dbh;
  gamma_slope_gi <- alpha*gamma_slope_dbh;
  gamma_aspect_gi <- alpha*gamma_aspect_dbh;
  gamma_twi_gi <- alpha*gamma_twi_dbh;
  for (i in 1:n_indiv) {
    alpha_indiv[i] <- beta_size_gi*tree_size_gi[i]
                      + beta_size_gi_squared*pow(tree_size_gi[i], 2.0)
                      + beta_basal_area_gi*basal_area_gi[i]
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
  beta_ppt_cool ~ normal(0.0, 1.0);
  beta_ppt_warm ~ normal(0.0, 1.0);
//  beta_tmn_cool ~ normal(0.0, 1.0);
//  beta_tmn_warm ~ normal(0.0, 1.0);
//  beta_TP_cool ~ normal(0.0, 1.0);
//  beta_TP_warm ~ normal(0.0, 1.0);
//  beta_pcool_sq ~ normal(0.0, 1.0);
//  beta_pwarm_sq ~ normal(0.0, 1.0);
//  beta_tcool_sq ~ normal(0.0, 1.0);
//  beta_twarm_sq ~ normal(0.0, 1.0);
//  beta_pcool_htr ~ normal(0.0, 1.0);
//  beta_tcool_htr ~ normal(0.0, 1.0);
//  beta_pwarm_htr ~ normal(0.0, 1.0);
//  beta_twarm_htr ~ normal(0.0, 1.0);
  beta_pcool_BA ~ normal(0.0, 1.0);
//  beta_tcool_BA ~ normal(0.0, 1.0);
  beta_pwarm_BA ~ normal(0.0, 1.0);
//  beta_twarm_BA ~ normal(0.0, 1.0);
//  beta_twarm_size ~ normal(0.0, 1.0);
  beta_pcool_size ~ normal(0.0, 1.0);
//  beta_tcool_size ~ normal(0.0, 1.0);
  beta_pwarm_size ~ normal(0.0, 1.0);
  sigma_indiv ~ normal(0.0, 1.0);
  sigma_site_gi ~ normal(0.0, 1.0);
  eta_sq ~ normal(0.0, 1.0);
  inv_rho_sq ~ normal(0.0, 1.0);
  sigma_sq ~ normal(0.0, 1.0);
  // priors for dbh component
  //
  beta_0_dbh ~ normal(0.0, 1.0);
  beta_size ~ normal(0.0, 1.0);
  beta_size_squared ~ normal(0.0, 1.0);
  beta_basal_area ~ normal(0.0, 1.0);
  beta_height_ratio ~ normal(0.0, 1.0);
  sigma_resid ~ normal(0.0, 1.0);
  sigma_site_dbh ~ normal(0.0, 1.0);
  sigma_species ~ normal(0.0, 1.0);
  // prior for shared components
  //
  // N.B.: alpha is defined with a lower bound of 0, so it is
  // effectively a half-normal prior
  //
  alpha ~ normal(0.0, sqrt(2.0));
  gamma_radiation_dbh ~ normal(0.0, 1.0);
  gamma_slope_dbh ~ normal(0.0, 1.0);
  gamma_aspect_dbh ~ normal(0.0, 1.0);
  gamma_twi_dbh ~ normal(0.0, 1.0);

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
  eps_species ~ normal(0.0, sigma_species);
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
