parameters {
  corr_matrix[2] omega;
}
transformed parameters {
  vector[2] zero;
  vector[2] tau;

  zero[1] <- 0.0;
  zero[2] <- 0.0;
  tau[1] <- 1.0;
  tau[2] <- 1.0;
}
model {
  omega ~ lkj_corr(2.0);
}
