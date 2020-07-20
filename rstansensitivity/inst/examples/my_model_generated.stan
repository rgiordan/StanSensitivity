data {
    real y;
  
  // Hyperparameters:
  real prior_mean;
  }
parameters {
  real mu;
  }
model {
  mu ~ normal(prior_mean, 1.0);
  y ~ normal(mu, 1.0);
  }

