data {
  int<lower=1> I;               // # items
  int<lower=1> J;               // # persons
  int<lower=1> N;               // # observations
  int<lower=1, upper=I> ii[N];  // item for n
  int<lower=1, upper=J> jj[N];  // person for n
  int<lower=0, upper=1> y[N];   // correctness for n
}
parameters {
  vector[J] theta;              // abilities
  vector[2] xi[I];              // alpha/beta pair vectors
  vector[2] mu;                 // vector for alpha/beta means
  vector<lower=0>[2] tau;       // vector for alpha/beta residual sds
  cholesky_factor_corr[2] L_Omega;

  // Hyperparameters:
  // Original values follow in comments.
  real lkj_concentration;  // 4
  real mu_loc; // 0
  real mu_1_scale; // 1
  real mu_2_scale; // 5
  real tau_loc; // 0.1
  real theta_loc; // 0
  real theta_scale; // 1
}
transformed parameters {
  vector[I] alpha;
  vector[I] beta;
  for (i in 1:I) {
    alpha[i] = exp(xi[i,1]);
    beta[i] = xi[i,2];
  }
}
model { target += 0; }generated quantities {
  corr_matrix[2] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}

