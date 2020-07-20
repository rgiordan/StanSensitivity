data {
    real y;
  }
parameters {
  real mu;
  
  // Hyperparameters:
  real prior_mean;
  }
model {
  mu ~ normal(prior_mean, 1.0);
  y ~ normal(mu, 1.0);
  }

