data {
  int<lower=0> N;
  vector[N] earn;
  int eth[N];
  vector[N] height;

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
model {
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] = a1[eth[i]] + a2[eth[i]] * height[i];

  // This is not really a prior on eta -- it is the likelihood specification
  // for a1 and a2 in the non-centered parameterization.
  eta1 ~ normal(0, 1);
  eta2 ~ normal(0, 1);

  mu_a1 ~ normal(mu_a1_loc, mu_a1_scale);
  mu_a2 ~ normal(mu_a2_loc, mu_a2_scale);

  // Tighten the random effect variance priors.
  sigma_a1 ~ normal(sigma_a1_loc, sigma_a1_scale);
  sigma_a2 ~ normal(sigma_a2_loc, sigma_a2_scale);
  sigma_y ~ cauchy(sigma_y_loc, sigma_y_scale);

  // To make the log posterior a density in log_sigma_*, add the log abs Jacobians
  // log(abs(d sigma_* / d log_sigma_*)) = log(abs(sigma_*)) = log_sigma_*
  target += log_sigma_y + log_sigma_a1 + log_sigma_a2;

  log_earn ~ normal(y_hat, sigma_y);
}

