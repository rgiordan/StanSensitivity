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
stan_sensitivity_list <- GetStanSensitivityModel(model_name, stan_data)
sens_time <- Sys.time()
sens_result <- GetStanSensitivityFromModelFit(sampling_result, draws_mat, stan_sensitivity_list)
sens_time <- Sys.time()- sens_time


##################################
# Inspect the results.

tidy_results <- GetTidyResult(draws_mat, sens_result)

dev.new()
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

dev.new()
ggplot(filter(tidy_results$sens_norm_df,
              abs(tidy_results$sens_norm_df$mean_sensitivity) > 1.0)) +
  geom_bar(aes(x=parameter, y=mean_sensitivity, fill=hyperparameter),
           stat="identity", position="dodge") +
  geom_errorbar(aes(x=parameter, ymin=lower_sensitivity,
                    ymax=upper_sensitivity, group=hyperparameter),
                position=position_dodge(0.9), width=0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


###########################
# Perturb and re-fit

perturb_par <- "R.3.3"
epsilon <- 0.1

# Importance sampling.
se_mean <- summary(sampling_result)$summary[, "se_mean"]
min_epsilon <- 2.0 * min(se_mean / abs(sens_result$sens_mat[perturb_par, ]))
if (epsilon < min_epsilon) {
  warning("The expected change is less than twice the mean standard error for every parameter.")
}

imp_list <- list()
num_eff_samples_list <- list()
eps_length <- 20
for (this_eps in seq(0, epsilon, length.out=eps_length)) {
  imp_sens_par_list <- stan_sensitivity_list$sens_par_list
  imp_sens_par_list[["R"]][3, 3] <- imp_sens_par_list[["R"]][3, 3] + this_eps

  imp_results <- GetImportanceSamplingFromModelFit(
    sampling_result, draws_mat, stan_sensitivity_list,
    imp_sens_par_list, lp_vec=sens_result$lp_vec)

  imp_diff <- imp_results$imp_means - colMeans(draws_mat)
  num_eff_samples_list[[length(num_eff_samples_list) + 1]] <-
    imp_results$eff_num_imp_samples
  imp_list[[length(imp_list) + 1]] <-
    data.frame(t(imp_diff)) %>%
    mutate(epsilon=this_eps) %>%
    melt(id.vars="epsilon") %>%
    rename(parameter=variable, imp_diff=value)
}


imp_df <- do.call(rbind, imp_list) %>%
  inner_join(filter(tidy_results$sens_df, hyperparameter=="R.3.3"), by="parameter")

max_eps <- max(unique(imp_df$epsilon)[unlist(num_eff_samples_list) > 500])

dev.new()
ggplot(filter(imp_df, parameter %in% sensitive_params[1:5])) +
  geom_point(aes(x=epsilon, y=imp_diff, color=parameter)) +
  geom_line(aes(x=epsilon, y=mean_sensitivity * epsilon, color=parameter)) +
  geom_vline(aes(xintercept=max_eps))

# Plot nonlinearities in the importance sampling.
dev.new()
ggplot(filter(imp_df, parameter %in% sensitive_params[21:30], grepl("alpha", parameter))) +
  geom_point(aes(x=epsilon, y=imp_diff, color=parameter)) +
  geom_line(aes(x=epsilon, y=mean_sensitivity * epsilon, color=parameter))


{
  cat("Imp. samp. difference:\t\t", imp_diff, "\n")
  cat("Sensitivity difference:\t\t", sens_result$sens_mat[perturb_par, , drop=FALSE] * epsilon, "\n")
  cat("Eff. # of importance samples:\t", imp_results$eff_num_imp_samples, "\n")
  # cat("Imp. samp. approx. err:\t\t", 2 * imp_se, "\n")
}


dev.new()
plot(imp_diff, sens_result$sens_mat[perturb_par, , drop=FALSE] * epsilon); abline(0, 1)

# Re-run MCMC
stan_data_perturb <- stan_data
stan_data_perturb[["R"]][3, 3] <- stan_data_perturb[["R"]][3, 3] + epsilon
mcmc_time <- Sys.time()
sampling_result_perturb <- sampling(model, data=stan_data_perturb, chains=1,
                                    iter=(num_samples + num_warmup_samples))
mcmc_time <- Sys.time() - mcmc_time
draws_mat_perturb <- extract(sampling_result_perturb, permute=FALSE)[,1,]
perturb_se <-
  data.frame(t(summary(sampling_result_perturb)$summary[, "se_mean"])) %>%
  mutate(method="mcmc") %>%
  melt(id.vars=c("method")) %>%
  rename(parameter=variable, se=value)

mcmc_diff <- colMeans(draws_mat_perturb) - colMeans(draws_mat)
mcmc_diff_df <-
  data.frame(t(mcmc_diff)) %>%
  mutate(method="mcmc") %>%
  melt(id.vars=c("method")) %>%
  rename(parameter=variable, mean_diff=value) %>%
  inner_join(filter(tidy_results$sens_df, hyperparameter==perturb_par), by="parameter") %>%
  inner_join(perturb_se, by=c("parameter", "method"))
head(mcmc_diff_df)

mcmc_diff_df <-
  mutate(mcmc_diff_df,
         top_change=abs(mean_diff) > quantile(abs(mean_diff), 0.9))

sensitive_params <-
  as.character(filter(tidy_results$sens_norm_df, abs(mean_sensitivity) > 1)$parameter)

dev.new()
ggplot(filter(mcmc_diff_df, parameter %in% sensitive_params)) +
  geom_point(aes(y=mean_diff, x=mean_sensitivity * epsilon, color=parameter), size=3) +
  geom_errorbar(aes(ymin=mean_diff - 2 * se, ymax=mean_diff + 2 * se, x=mean_sensitivity * epsilon)) +
  geom_errorbarh(aes(y=mean_diff, x=mean_sensitivity * epsilon,
                     xmin=lower_sensitivity * epsilon, xmax=upper_sensitivity * epsilon)) +
  geom_abline(aes(slope=1, intercept=0))

# Omega.2.2 is off.
top_misses <-
  arrange(mcmc_diff_df, desc(abs(mean_diff - epsilon * mean_sensitivity)))[1:10, "parameter"]

dev.new()
ggplot(filter(mcmc_diff_df, parameter %in% top_misses)) +
  geom_point(aes(y=mean_diff, x=mean_sensitivity * epsilon, color=parameter), size=3) +
  geom_errorbar(aes(ymin=mean_diff - 2 * se, ymax=mean_diff + 2 * se, x=mean_sensitivity * epsilon)) +
  geom_errorbarh(aes(y=mean_diff, x=mean_sensitivity * epsilon,
                     xmin=lower_sensitivity * epsilon, xmax=upper_sensitivity * epsilon)) +
  geom_abline(aes(slope=1, intercept=0))

{
  cat("MCMC difference:\t\t", colMeans(draws_mat_perturb) - colMeans(draws_mat), "\n")
  cat("Imp. samp. difference:\t\t", imp_diff, "\n")
  cat("Sensitivity difference:\t\t", sens_mat[perturb_par, , drop=FALSE] * epsilon, "\n")
  cat("Eff. # of importance samples:\t", imp_results$eff_num_imp_samples, "\n")
  cat("Imp. samp. approx. err:\t\t", 2 * imp_se, "\n")
  cat("MCMC approx. err:\t\t", 2 * perturb_se, "\n")
}
