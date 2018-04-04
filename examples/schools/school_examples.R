library(rstan)
library(rstansensitivity)

library(ggplot2)
library(dplyr)
library(reshape2)

rstan_options(auto_write=TRUE)

# Set this to be the appropriate location of the repo on your computer.
# Run from anywhere in the StanSensitivity repository.
git_repo <- system("git rev-parse --show-toplevel", intern=TRUE)
example_dir <- file.path(git_repo, "examples/schools/")
model_name <- GenerateSensitivityFromModel(
  file.path(example_dir, "models/schools-1.stan"))

##################################
# Compile and run the base model.

model <- stan_model(GetSamplingModelFilename(model_name))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(example_dir, "schools-1.data.R", sep=""), local=stan_data)
stan_data <- as.list(stan_data)

# For now, you must use chains=1 for now to avoid confusion around get_inits.
num_warmup_samples <- 2500
num_samples <- 2500
sampling_file <- paste(model_name, "_sampling.Rdata", sep="")
if (!file.exists(sampling_file)) {
  print("Running sampler.")
  sampling_time <- time.time()
  sampling_result <- sampling(model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples))
  sampling_time <- time.time() - sampling_time
  save(sampling_result, sampling_result, sampling_time,
       file=sampling_file)
} else {
  print(sprintf("Loading cached samples from %s", sampling_file))
  load(sampling_file)  
}
print(summary(sampling_result))

stan_sensitivity_list <- GetStanSensitivityModel(model_name, stan_data)
sens_result <- GetStanSensitivityFromModelFit(sampling_result, stan_sensitivity_list)

tidy_results <- GetTidyResult(sens_result)
PlotSensitivities(filter(tidy_results, abs(normalized_sensitivity)  > 1.0))
PlotSensitivities(filter(tidy_results, grepl("beta", parameter)))


###########################
# Perturb and re-fit to check the sensitivity measurements.

perturb_par <- "R.3.3"
epsilon <- 0.1

# Check that the perturbation is big enough that we expect to produce
# a measurable difference in the output.
se_mean <- summary(sampling_result)$summary[, "se_mean"]
min_epsilon <- 2.0 * min(se_mean / abs(sens_result$sens_mat[perturb_par, ]))
if (epsilon < min_epsilon) {
  warning("The expected change is less than twice the mean standard error for every parameter.")
}


# Re-run MCMC
stan_data_perturb <- stan_data
stan_data_perturb[["R"]][3, 3] <- stan_data_perturb[["R"]][3, 3] + epsilon
perturbed_sampling_file <- paste(model_name, "_perturbed_sampling.Rdata", sep="")
if (!file.exists(perturbed_sampling_file)) {
  print("Running sampler.")
  sampling_time_perturb <- time.time()
  sampling_result_perturb <- sampling(model, data=stan_data_perturb, chains=1,
                                      iter=(num_samples + num_warmup_samples))
  sampling_time_perturb <- time.time() - sampling_time_perturb
  save(sampling_result_perturb, stan_data_perturb, sampling_time_perturb,
       file=perturbed_sampling_file)
} else {
  print(sprintf("Loading cached perturbed samples from %s", perturbed_sampling_file))
  load(perturbed_sampling_file)  
}

draws_mat_perturb <- extract(sampling_result_perturb, permute=FALSE)[,1,]

perturb_se <-
  data.frame(t(summary(sampling_result_perturb)$summary[, "se_mean"])) %>%
  mutate(method="mcmc") %>%
  melt(id.vars=c("method")) %>%
  rename(parameter=variable, se=value)

mcmc_diff <- colMeans(draws_mat_perturb) - colMeans(sens_result$draws_mat)
mcmc_diff_df <-
  data.frame(t(mcmc_diff)) %>%
  mutate(method="mcmc") %>%
  melt(id.vars=c("method")) %>%
  rename(parameter=variable, mean_diff=value) %>%
  inner_join(filter(tidy_results, hyperparameter==perturb_par), by="parameter") %>%
  inner_join(perturb_se, by=c("parameter", "method"))
head(mcmc_diff_df)

mcmc_diff_df <-
  mutate(mcmc_diff_df,
         top_change=abs(mean_diff) > quantile(abs(mean_diff), 0.9))

sensitive_params <-
  as.character(filter(tidy_results, abs(sensitivity) > 1)$parameter)

# Compare the actual effects to the predicted effects based on a linear approximation.
ggplot(filter(mcmc_diff_df, parameter %in% sensitive_params)) +
  geom_point(aes(y=mean_diff, x=sensitivity * epsilon, color=parameter), size=3) +
  geom_errorbar(aes(ymin=mean_diff - 2 * se, ymax=mean_diff + 2 * se, x=sensitivity * epsilon)) +
  geom_errorbarh(aes(y=mean_diff, x=sensitivity * epsilon,
                     xmin=(sensitivity - 2 * sensitivity_se) * epsilon,
                     xmax=(sensitivity + 2 * sensitivity_se) * epsilon)) +
  geom_abline(aes(slope=1, intercept=0))
