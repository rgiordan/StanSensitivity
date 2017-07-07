
#########################################
# Importance sampling.

epsilon <- 0.01

imp_sens_par_list <- stan_sensitivity_list$sens_par_list
imp_sens_par_list$y_var <- imp_sens_par_list$y_var + epsilon
imp_results <- GetImportanceSamplingFromModelFit(
  sampling_result, draws_mat, stan_sensitivity_list,
  imp_sens_par_list, lp_vec=lp_vec)

cat("Effective number of importance samples: ", imp_results$eff_num_imp_samples, "\n")

if (FALSE) {
  hist(log10(imp_results$imp_weights), 100)
  abline(v=-log10(nrow(draws_mat)))
}

imp_mu_diff <- imp_results$imp_means["mu"] - colMeans(draws_mat[, "mu", drop=FALSE])

eff_num_samples <- summary(sampling_result)$summary[, "n_eff"] 
imp_se <- summary(sampling_result)$summary[, "se_mean"] *
          sqrt(eff_num_samples / imp_results$eff_num_imp_samples)

cat("Importance sampling difference:\t\t", imp_mu_diff, "\n")
cat("Sensitivity difference:\t\t\t", sens_mat["y_var", "mu"] * epsilon, "\n")



#########################################
# Perturb and re-draw.

stan_data_perturb <- stan_data
stan_data_perturb$y_var <- stan_data_perturb$y_var + epsilon
sampling_result_perturb <- sampling(model, data=stan_data_perturb, chains=1,
                                    iter=(num_samples + num_warmup_samples))
draws_mat_perturb <- extract(sampling_result_perturb, permute=FALSE)[,1,]

print("MCMC difference: ")
print(mean(draws_mat_perturb[, "mu"]) - mean(draws_mat[, "mu"]))
print("Importance sampling difference: ")
print(imp_mu_diff)
print("Sensitivity difference: ")
sens_mat["y_var", "mu"] * epsilon


