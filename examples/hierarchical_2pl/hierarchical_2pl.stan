data {
  int<lower=1> I;               // # items
  int<lower=1> J;               // # persons
  int<lower=1> N;               // # observations
  int<lower=1, upper=I> ii[N];  // item for n
  int<lower=1, upper=J> jj[N];  // person for n
  int<lower=0, upper=1> y[N];   // correctness for n
}
hyperparameters {
  // Original values follow in comments.
  real lkj_concentration;  // 4
  real mu_loc; // 0
  real mu_1_scale; // 1
  real mu_2_scale; // 5
  real tau_loc; // 0.1
  real theta_loc; // 0
  real theta_scale; // 1
}
parameters {
  vector[J] theta;              // abilities
  vector[2] xi[I];              // alpha/beta pair vectors
  vector[2] mu;                 // vector for alpha/beta means
  vector<lower=0>[2] tau;       // vector for alpha/beta residual sds
  cholesky_factor_corr[2] L_Omega;
}
transformed parameters {
  vector[I] alpha;
  vector[I] beta;
  for (i in 1:I) {
    alpha[i] = exp(xi[i,1]);
    beta[i] = xi[i,2];
  }
}
model {
  matrix[2,2] L_Sigma;
  L_Sigma = diag_pre_multiply(tau, L_Omega);
  for (i in 1:I)
    xi[i] ~ multi_normal_cholesky(mu, L_Sigma);
  theta ~ normal(theta_loc, theta_scale);
  L_Omega ~ lkj_corr_cholesky(lkj_concentration);
  mu[1] ~ normal(mu_loc, mu_1_scale);
  tau[1] ~ exponential(tau_loc);
  mu[2] ~ normal(mu_loc, mu_2_scale);
  tau[2] ~ exponential(tau_loc);
  y ~ bernoulli_logit(alpha[ii] .* (theta[jj] - beta[ii]));
}
generated quantities {
  corr_matrix[2] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}
