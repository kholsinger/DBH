data {
  int<lower=0> n_obs;
  int<lower=0> n_plots;
  int<lower=0> n_species;
  vector[n_obs] dbh_2;
  vector[n_obs] dbh_1;
  vector[n_obs] tree_size;
  vector[n_obs] height_ratio;
  vector[n_obs] radiation;
  vector[n_obs] slope;
  vector[n_obs] aspect;
  vector[n_obs] twi;
  int<lower=0> plot[n_obs];
  int<lower=0> species[n_obs];
}
parameters {
  real beta_0;
  real beta_size;
  real beta_height_ratio;
  real gamma_radiation;
  real gamma_slope;
  real gamma_aspect;
  real gamma_twi;
  real<lower=0> sigma_resid;
  real<lower=0> sigma_plot;
  real<lower=0> sigma_species;
  vector[n_plots] eps_plot;
  vector[n_species] eps_species;
}
transformed parameters {
  vector[n_obs] mu_indiv;
  vector[n_obs] log_mu_dbh_inc;

  for (i in 1:n_obs) {
    log_mu_dbh_inc[i] <- beta_size*tree_size[i]
                         + beta_height_ratio*height_ratio[i]
                         + gamma_radiation*radiation[plot[i]]
                         + gamma_slope*slope[plot[i]]
                         + gamma_aspect*aspect[plot[i]]
                         + gamma_twi*twi[plot[i]]
                         + eps_species[species[i]]
                         + eps_plot[plot[i]];
  }
  mu_indiv <- dbh_1 + exp(log_mu_dbh_inc);
}
model {
  beta_0 ~ normal(0.0, 1.0);
  beta_size ~ normal(0.0, 1.0);
  beta_height_ratio ~ normal(0.0, 1.0);
  gamma_radiation ~ normal(0.0, 1.0);
  gamma_slope ~ normal(0.0, 1.0);
  gamma_aspect ~ normal(0.0, 1.0);
  gamma_twi ~ normal(0.0, 1.0);
  sigma_resid ~ normal(0.0, 1.0);
  sigma_plot ~ normal(0.0, 1.0);
  sigma_species ~ normal(0.0, 1.0);

  eps_species ~ normal(beta_0, sigma_species);
  eps_plot ~ normal(0.0, sigma_plot);
  dbh_2 ~ normal(mu_indiv, sigma_resid);
}