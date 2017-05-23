data {
real U;
int<lower=0> N_censored;
int<lower=0> N_observed;
real<upper=U> y[N_observed];
}
parameters {
real mu;real weights[N_observed];
real y_var;

}
model {
for (n in 1:N_observed) {
  target += weights[n] * (
      normal_lpdf(y[n] | mu, y_var) - normal_lcdf(U | mu, y_var));
  }
target += N_censored * log1m(normal_cdf(U, mu, y_var));
}

