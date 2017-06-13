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
}
parameters {
  real beta[8];
  real<lower=0> theta;
  real phi;
  vector[3] alpha[M];
  vector[3] gamma;
  cov_matrix[3] Omega;

  // Hyperparameters:
  matrix[3, 3] R;
  real wishart_df;
  real hp_mean;
  real hp_var;
}
model {
  target += 0.;
}
