data {
  int<lower=1> N;
  int<lower=0> y[N];

  // Hyperparameters:
  real cauchy_loc_alpha;
  real cauchy_loc_beta;
  real cauchy_scale_alpha;
  real cauchy_scale_beta;
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
}
transformed parameters {
  real log_alpha;
  real log_beta;
  log_alpha = log(alpha);
  log_beta = log(beta);
}
model {
  alpha ~ cauchy(cauchy_loc_alpha, cauchy_scale_alpha);
  beta ~ cauchy(cauchy_loc_beta, cauchy_scale_beta);
  for (i in 1:N)
    target += neg_binomial_lpmf(y[i] | alpha, beta);
}

