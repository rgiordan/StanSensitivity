# Compute the variance of the covariances and their ratios.

library(devtools)
# install_local("/home/rgiordan/Documents/git_repos/StanSensitivity/rstansensitivity", force=TRUE)

boot_results <- BootstrapSensitivityMatrix(draws_mat, sens_result$grad_mat, num_boot=100)

sens_mat_quantiles <-
  GetArrayQuantiles(boot_results$sens_mat_array, alpha=0.1)
sens_mat_normalized_quantiles <-
  GetArrayQuantiles(boot_results$sens_mat_normalized_array, alpha=0.1)


param_names <- stan_sensitivity_list$param_names
sens_param_names <- stan_sensitivity_list$sens_param_names
grad_mat <- sens_result$grad_mat

rownames(grad_mat) <- sens_param_names
grad_mat <- grad_mat[setdiff(sens_param_names, param_names), , drop=FALSE]

sens_mat <- GetSensitivityFromGrads(grad_mat, draws_mat)


sens_mat_quantiles$lower["y_var", ]
sens_mat_quantiles$upper["y_var", ]
sens_mat["y_var", ]

