if (FALSE) {
  library(devtools)
  install_local("/home/rgiordan/Documents/git_repos/StanSensitivity/rstansensitivity", force=TRUE)  
}

#########################################
# Importance sampling to check linearity.


imp_sens_par_list <- stan_sensitivity_list$sens_par_list
if (grepl("normal_censored.stan$", base_model_name)) {
  perturb_par <- "y_var"
  epsilon <- 0.03
} else if (grepl("negative_binomial.stan$", base_model_name)) {
  epsilon <- 2.0
  perturb_par <- "cauchy_scale_alpha"
}

param_names <- stan_sensitivity_list$param_names
se_mean <- summary(sampling_result)$summary[param_names, "se_mean"]
min_epsilon <- 2.0 * min(se_mean / abs(sens_mat[perturb_par, param_names]))
if (epsilon < min_epsilon) {
  warning("The expected change is less than twice the mean standard error for every parameter.")
}
imp_sens_par_list[[perturb_par]] <- imp_sens_par_list[[perturb_par]] + epsilon

imp_results <- GetImportanceSamplingFromModelFit(
  sampling_result, draws_mat, stan_sensitivity_list,
  imp_sens_par_list, lp_vec=sens_result$lp_vec)


if (FALSE) {
  hist(log10(imp_results$imp_weights), 100)
  abline(v=-log10(nrow(draws_mat)), lwd=3, col="red")
}

imp_diff <- imp_results$imp_means - colMeans(draws_mat)

eff_num_samples <- summary(sampling_result)$summary[, "n_eff"] 
imp_se <- summary(sampling_result)$summary[, "se_mean"] *
          sqrt(nrow(draws_mat) / imp_results$eff_num_imp_samples)

summary(sampling_result)
{
  cat("Effective number of importance samples: ", imp_results$eff_num_imp_samples, "\n")
  cat("Importance sampling difference:\t\t", imp_diff, "\n")
  cat("Sensitivity difference:\t\t\t", sens_mat[perturb_par, , drop=FALSE] * epsilon, "\n")
}



#########################################
# Perturb and re-draw.

stan_data_perturb <- stan_data
stan_data_perturb[[perturb_par]] <- stan_data_perturb[[perturb_par]] + epsilon
sampling_result_perturb <- sampling(model, data=stan_data_perturb, chains=1,
                                    iter=(num_samples + num_warmup_samples))
draws_mat_perturb <- extract(sampling_result_perturb, permute=FALSE)[,1,]
perturb_se <- summary(sampling_result_perturb)$summary[, "se_mean"]

{
  cat("MCMC difference:\t\t", colMeans(draws_mat_perturb) - colMeans(draws_mat), "\n")
  cat("Imp. samp. difference:\t\t", imp_diff, "\n")
  cat("Sensitivity difference:\t\t", sens_mat[perturb_par, , drop=FALSE] * epsilon, "\n")
  cat("Eff. # of importance samples:\t", imp_results$eff_num_imp_samples, "\n")
  cat("Imp. samp. approx. err:\t\t", 2 * imp_se, "\n")
  cat("MCMC approx. err:\t\t", 2 * perturb_se, "\n")
}

eps_max <- 0.2
eps_len <- 5
eps_range <- seq(0.0, eps_max, length.out=eps_len)
mcmc_diff_mat <- matrix(NA, length(eps_range), ncol(draws_mat))
for (n in 1:length(eps_range)) {
  stan_data_perturb <- stan_data
  stan_data_perturb[[perturb_par]] <- stan_data_perturb[[perturb_par]] + eps_range[n]
  sampling_result_perturb <- sampling(model, data=stan_data_perturb, chains=1,
                                      iter=(num_samples + num_warmup_samples))
  draws_mat_perturb <- extract(sampling_result_perturb, permute=FALSE)[,1,]
  mcmc_diff_mat[n, ] <- colMeans(draws_mat_perturb) - colMeans(draws_mat)
}

sens_df <-
  SensitivityMatrixToDataframe(sens_mat, rownames(sens_mat), colnames(sens_mat)) %>%
  filter(parameter != "lp__", hyperparameter == perturb_par)

colnames(mcmc_diff_mat) <- colnames(draws_mat)
mcmc_diff_df <-
  data.frame(mcmc_diff_mat) %>%
  mutate(epsilon=eps_range) %>%
  melt(id.var="epsilon") %>%
  filter(variable != "lp__") %>%
  mutate(method="mcmc")

ggplot() +
  geom_point(data=mcmc_diff_df, aes(x=epsilon, y=value, color=variable)) +
  geom_abline(data=sens_df, aes(intercept=0, slope=value, color=parameter), lwd=2) +
  expand_limits(x=0, y=0)

ggplot(sens_df) +
  geom_abline(aes(intercept=0, slope=value, color=parameter, linetype=hyperparameter))
