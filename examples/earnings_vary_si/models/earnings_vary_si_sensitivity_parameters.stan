data {
  int<lower=0> N;
  vector[N] earn;
  int eth[N];
  vector[N] height;
}
transformed data {
  vector[N] log_earn;

  log_earn = log(earn);
}
parameters {
  vector[4] a1;
  vector[4] a2;
  real<lower=0> sigma_a1;
  real<lower=0> sigma_a2;
  real<lower=0> sigma_y;
  real mu_a1;
  real mu_a2;

  // Hyperparameters:
  real mu_a1_loc;
  real mu_a1_scale;
  real mu_a2_loc;
  real mu_a2_scale;
  real sigma_a1_loc;
  real sigma_a1_scale;
  real sigma_a2_loc;
  real sigma_a2_scale;
  real sigma_y_loc;
  real sigma_y_scale;
}
model { target += 0; }
