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
model {
  vector[N] y_hat;

  for (i in 1:N)
    y_hat[i] = a1[eth[i]] + a2[eth[i]] * height[i];

  mu_a1 ~ normal(mu_a1_loc, mu_a1_scale);
  mu_a2 ~ normal(mu_a2_loc, mu_a2_scale);
  a1 ~ normal(10 * mu_a1, sigma_a1);
  a2 ~ normal(0.01 * mu_a2, sigma_a2);
  sigma_a1 ~ cauchy(sigma_a1_loc, sigma_a1_scale);
  sigma_a2 ~ cauchy(sigma_a2_loc, sigma_a2_scale);
  sigma_y ~ cauchy(sigma_y_loc, sigma_y_scale);
  log_earn ~ normal(y_hat, sigma_y);
}

