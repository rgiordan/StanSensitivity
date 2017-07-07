data {
  int<lower=1> N;
  int<lower=0> y[N];
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;

  // Hyperparameters:
  real weights[N];
  real cauchy_loc_alpha;
  real cauchy_loc_beta;
  real cauchy_scale_alpha;
  real cauchy_scale_beta;
}
model { target += 0; }
