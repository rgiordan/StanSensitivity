// This model has a closed form mean as a function of the weights and
// so can be used as a unit test for derivatives.
data {
  int<lower=1> n_obs;
  vector[n_obs] x;
  vector[n_obs] w;
  real<lower=0> x_sd;
  real<lower=0> mu_sd;
  real mu_mean;
}
parameters {
  real mu;
}
transformed parameters {
  vector[n_obs] log_lik;
  for (n in 1:n_obs) {
    log_lik[n] = normal_lpdf(x[n] | mu, x_sd);
  }
}
model {
  mu ~ normal(mu_mean, mu_sd);
  target += dot_product(w, log_lik);
}
