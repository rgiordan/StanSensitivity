# This is intended to be run after run_examples.

GetPosteriorMeans <- function(stan_data) {
  result <- sampling(model, data=stan_data, chains=1, iter=num_samples * 2)
  return(summary(result)$summary[,"mean"])
}


epsilon <- 1e-3
weight_rows <- grepl("weights", sens_param_names)
weight_param_names <- sens_param_names[weight_rows]
param_means <- summary(result)$summary[, "mean"]

pert_mat <- matrix(NA, sum(weight_rows), ncol(sens_mat))
rownames(pert_mat) <- weight_param_names
colnames(pert_mat) <- colnames(sens_mat)

for (wi in 1:length(weight_param_names)) {
  stan_data_perturb <- stan_data
  stan_data_perturb$weights[wi] <- stan_data_perturb$weights[wi] + epsilon
  pert_mat[wi, ] <- GetPosteriorMeans(stan_data_perturb) - param_means
}

pert_mat
sens_mat[weight_param_names, ] * epsilon
summary(result)$summary
