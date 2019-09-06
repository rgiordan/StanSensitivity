// Test a model with multiple parameters for reweighting.
data {
  int<lower=1> n_obs;
  int<lower=1> x_dim;
  matrix[n_obs, x_dim] x;
  vector[n_obs] w;
}
parameters {
  vector[x_dim] mu;
  real<lower=0> sigma;
}
transformed parameters {
  vector[n_obs] log_lik;
  real log_sigma;
  log_sigma = log(sigma);
  for (n in 1:n_obs) {
    log_lik[n] = normal_lpdf(x[n, 1:x_dim] | mu, sigma);
  }
}
model {
  target += dot_product(w, log_lik);
}
