data {
  int<lower=1> N;
  int<lower=0> y[N];
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;

  // Hyperparameters:
  real cauchy_loc_alpha;
  real cauchy_loc_beta;
  real cauchy_scale_alpha;
  real cauchy_scale_beta;
}
transformed parameters {
  real log_alpha;
  real log_beta;
  log_alpha = log(alpha);
  log_beta = log(beta);
}
model { target += 0; }
