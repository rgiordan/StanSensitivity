data {
  int<lower=0> N;
  int<lower=0> M;
  vector[N] LRT;
  int school[N];
  int School_denom[N, 3];
  int School_gender[N, 2];
  int VR[N, 2];
  real Y[N];
  int Gender[N];

  // Hyperparameters:
  matrix[3, 3] R;
  real wishart_df;
  real hp_mean;
  real hp_var;
}
parameters {
  real beta[8];
  real<lower=0> theta;
  real phi;
  vector[3] alpha[M];
  vector[3] gamma;
  cov_matrix[3] Omega;
}
model {
  real Ymu[N];
  for(p in 1:N)
    Ymu[p] = alpha[school[p],1]
      + alpha[school[p],2] * LRT[p]
      + alpha[school[p],3] * VR[p,1]
      + beta[1] * LRT[p] * LRT[p]
      + beta[2] * VR[p,2]
      + beta[3] * Gender[p]
      + beta[4] * School_gender[p,1]
      + beta[5] * School_gender[p,2]
      + beta[6] * School_denom[p,1]
      + beta[7] * School_denom[p,2]
      + beta[8] * School_denom[p,3];

  Y ~ normal(Ymu,  exp(-0.5 * (theta + phi * LRT)));

  # Priors for fixed effects:
  beta ~ normal(hp_mean, hp_var);
  theta ~ normal(hp_mean, hp_var);
  phi ~ normal(hp_mean, hp_var);

  # Priors for random coefficients:
  alpha ~ multi_normal_prec(gamma, Omega);

  # Hyper-priors:
  gamma ~ normal(hp_mean, hp_var);
  Omega ~ wishart(wishart_df, R);
}

