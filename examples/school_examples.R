library(rstan)
library(rstansensitivity)

library(ggplot2)
library(dplyr)
library(reshape2)

rstan_options(auto_write=TRUE)

# Set this to be the appropriate location of the repo on your computer.
example_directory <- file.path(
  Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/examples/example_models")

base_model_name <- file.path(example_directory, "schools/schools-1.stan")
python_script <- file.path(Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/python/generate_models.py")
model_name <- GenerateSensitivityFromModel(base_model_name, python_script=python_script)


##################################
# Compile and run the base model.

model <- stan_model(paste(model_name, "_generated.stan", sep=""))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(model_name, "data.R", sep="."), local=stan_data)
stan_data <- as.list(stan_data)

# For now, you must use chains=1 for now to avoid confusion around get_inits.
# The script currently assumes the same number of warm-up draws as final samples.
num_warmup_samples <- 1000
num_samples <- 1000
sampling_result <- sampling(model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples))
print(summary(sampling_result))

##################################
# Get the sensitivity model and sensitivity.

draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
stan_sensitivity_list <- GetStanSensitivityModel(sampling_result, model_name, stan_data)
sens_time <- Sys.time()
sens_result <- GetStanSensitivityFromModelFit(sampling_result, draws_mat, stan_sensitivity_list)
sens_time <- Sys.time()- sens_time


##################################
# Inspect the results.

tidy_results <- GetTidyResult(draws_mat, sens_result)  

ggplot(filter(tidy_results$sens_norm_df,
              abs(mean_sensitivity) > 1.0)) +
  geom_bar(aes(x=parameter, y=mean_sensitivity, fill=hyperparameter),
           stat="identity", position="dodge") +
  geom_errorbar(aes(x=parameter, ymin=lower_sensitivity,
                    ymax=upper_sensitivity, group=hyperparameter),
                position=position_dodge(0.9), width=0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(filter(tidy_results$sens_norm_df,
              grepl("beta", parameter))) +
  geom_bar(aes(x=parameter, y=mean_sensitivity, fill=hyperparameter),
           stat="identity", position="dodge") +
  geom_errorbar(aes(x=parameter, ymin=lower_sensitivity,
                    ymax=upper_sensitivity, group=hyperparameter),
                position=position_dodge(0.9), width=0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(filter(tidy_results$sens_df, !grepl("weight", hyperparameter))) +
  geom_bar(aes(x=parameter, y=mean_sensitivity, fill=hyperparameter),
           stat="identity", position="dodge") +
  geom_errorbar(aes(x=parameter, ymin=lower_sensitivity,
                    ymax=upper_sensitivity, group=hyperparameter),
                position=position_dodge(0.9), width=0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


###########################
# Perturb and re-fit

perturb_par <- "R.3.3"
epsilon <- 0.2

# Importance sampling.
# Note: it appears that param names doesn't match up with the summary names
# if you have covariance matrices.
# param_names <- stan_sensitivity_list$param_names
# se_mean <- summary(sampling_result)$summary[param_names, "se_mean"]
# min_epsilon <- 2.0 * min(se_mean / abs(sens_mat[perturb_par, param_names]))
# if (epsilon < min_epsilon) {
#   warning("The expected change is less than twice the mean standard error for every parameter.")
# }
imp_sens_par_list <- stan_sensitivity_list$sens_par_list
imp_sens_par_list[["R"]][3, 3] <- imp_sens_par_list[["R"]][3, 3] + epsilon

imp_results <- GetImportanceSamplingFromModelFit(
  sampling_result, draws_mat, stan_sensitivity_list,
  imp_sens_par_list, lp_vec=sens_result$lp_vec)

imp_diff <- imp_results$imp_means - colMeans(draws_mat)

{
  cat("Imp. samp. difference:\t\t", imp_diff, "\n")
  cat("Sensitivity difference:\t\t", sens_result$sens_mat[perturb_par, , drop=FALSE] * epsilon, "\n")
  cat("Eff. # of importance samples:\t", imp_results$eff_num_imp_samples, "\n")
  # cat("Imp. samp. approx. err:\t\t", 2 * imp_se, "\n")
}

plot(imp_diff, sens_result$sens_mat[perturb_par, , drop=FALSE] * epsilon)

# Re-run MCMC
stan_data_perturb <- stan_data
stan_data_perturb[["R"]][3, 3] <- stan_data_perturb[["R"]][3, 3] + epsilon
sampling_result_perturb <- sampling(model, data=stan_data_perturb, chains=1,
                                    iter=(num_samples + num_warmup_samples))
draws_mat_perturb <- extract(sampling_result_perturb, permute=FALSE)[,1,]
perturb_se <- summary(sampling_result_perturb)$summary[, "se_mean"]
mcmc_diff <- colMeans(draws_mat_perturb) - colMeans(draws_mat)

{
  cat("MCMC difference:\t\t", colMeans(draws_mat_perturb) - colMeans(draws_mat), "\n")
  cat("Imp. samp. difference:\t\t", imp_diff, "\n")
  cat("Sensitivity difference:\t\t", sens_mat[perturb_par, , drop=FALSE] * epsilon, "\n")
  cat("Eff. # of importance samples:\t", imp_results$eff_num_imp_samples, "\n")
  cat("Imp. samp. approx. err:\t\t", 2 * imp_se, "\n")
  cat("MCMC approx. err:\t\t", 2 * perturb_se, "\n")
}


plot(mcmc_diff, sens_result$sens_mat[perturb_par, , drop=FALSE] * epsilon)



