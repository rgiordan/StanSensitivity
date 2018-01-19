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
  vector[4] eta1;
  vector[4] eta2;
  real mu_a1;
  real mu_a2;
  real log_sigma_a1;
  real log_sigma_a2;
  real log_sigma_y;

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
transformed parameters {
  real<lower=0> sigma_a1;
  real<lower=0> sigma_a2;
  real<lower=0> sigma_y;

  vector[4] a1;
  vector[4] a2;

  sigma_a1 = exp(log_sigma_a1);
  sigma_a2 = exp(log_sigma_a2);
  sigma_y = exp(log_sigma_y);

  a1 = 10 * mu_a1 + sigma_a1 * eta1;
  a2 = 0.1 * mu_a2 + sigma_a2 * eta2;
}
model { target += 0; }
