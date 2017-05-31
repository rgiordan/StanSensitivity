# This is intended to be run after run_examples.

GetPosteriorMeans <- function(stan_data) {
  result <- sampling(model, data=stan_data, chains=1, iter=num_samples * 2)
  return(summary(result)$summary[,"mean"])
}

hyperparam_rows <- !grepl("weights", rownames(sens_mat))
hyperparam_names <- rownames(sens_mat)[param_rows]
param_means <- summary(result)$summary[, "mean"]

pert_mat <- matrix(NA, length(hyperparam_names), ncol(sens_mat))
rownames(pert_mat) <- hyperparam_names
colnames(pert_mat) <- colnames(sens_mat)
param_names <- setdiff(colnames(pert_mat), "lp__")

# Set the perturbation well outside the standard error range.
epsilon <- 10 * min(summary(result)$summary[param_names, "se_mean"])

for (i in 1:length(hyperparam_names)) {
  hp_name <- hyperparam_names[i]
  stan_data_perturb <- stan_data
  stan_data_perturb[[hp_name]] <- stan_data_perturb[[hp_name]]  + epsilon
  pert_mat[i, ] <- GetPosteriorMeans(stan_data_perturb) - param_means
}

pert_mat[, param_names]
sens_mat[hyperparam_names, param_names] * epsilon
