data {
  real U;
  int<lower=0> N_censored;
  int<lower=0> N_observed;
  real<upper=U> y[N_observed];
}
parameters {
  real mu;

  // Hyperparameters:
  real weights[N_observed];
  real y_var;
}
model { target += 0; }
