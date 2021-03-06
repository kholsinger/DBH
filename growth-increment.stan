data {
  int<lower=0> n_obs;
  int<lower=0> n_years;
  int<lower=0> n_sites;
  int<lower=0> n_indiv;
  int<lower=0> n_months;
  matrix[n_years,n_months] ppt;
  matrix[n_years,n_months] tmn;
  vector[n_obs] gi;
  int<lower=0> year[n_obs];
  int<lower=0> indiv[n_obs];
  int<lower=0> site[n_indiv];
  // for index calculations
  vector[n_months] ppt_mean;
  vector[n_months] tmn_mean;
  // for robust regression
  real nu;
}
parameters {
  vector[n_months] beta_ppt;
  vector[n_months] beta_tmn;
  vector[n_indiv] mu_indiv;
  real<lower=0> sigma_resid;
  real<lower=0> sigma_indiv;
}
transformed parameters {
  matrix[n_years,n_indiv] mu_year_indiv;
  vector[n_years] mu_year;

  // n_months is number of months of prior weather included in
  // calculating expectation
  //
  mu_year <- ppt*beta_ppt + tmn*beta_tmn;
  // expectation for individual j in year i is sum of year
  // and indivdidual effects
  //
  for (i in 1:n_years) {
    for (j in 1:n_indiv) {
      mu_year_indiv[i,j] <- mu_year[i] + mu_indiv[j];
    }
  }
}
model {
  // priors
  //
  beta_ppt ~ normal(0.0, 1.0);
  beta_tmn ~ normal(0.0, 1.0);
  sigma_resid ~ cauchy(0.0, 2.5);
  sigma_indiv ~ cauchy(0.0, 2.5);

  // likelihood
  //
  for (i in 1:n_indiv) {
    mu_indiv[i] ~ normal(0.0, sigma_indiv);
  }
  // individual site x year combinations
  //
  for (i in 1:n_obs) {
    gi[i] ~ student_t(nu, mu_year_indiv[year[i],indiv[i]], sigma_resid);
  }
}
generated quantities {
  real idx_site;
  vector[n_obs] log_lik;

  idx_site <- ppt_mean'*beta_ppt + tmn_mean'*beta_tmn;
  // calculate log likelihoods
  //
  for (i in 1:n_obs) {
    log_lik[i] <- student_t_log(gi[i], nu,
                                mu_year_indiv[year[i],indiv[i]],
                                sigma_resid);
  }
}
