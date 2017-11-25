library(rstan)
library(rstansensitivity)

rstan_options(auto_write=TRUE)

# Set this to be the appropriate location of the repo on your computer.
example_directory <- file.path(
  Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/examples/example_models")

if (FALSE) {
  base_model_name <- file.path(
      example_directory, "negative_binomial/negative_binomial.stan")
  num_warmup_samples <- 5000
  num_samples <- 5000
}


if (TRUE) {
  base_model_name <- file.path(
      example_directory, "normal_censored/normal_censored.stan")
  num_warmup_samples <- 5000
  num_samples <- 5000
}

# Set this to be the location of the python script, or run it by hand following
# the directions in the README.
python_script <- file.path(Sys.getenv("GIT_REPO_LOC"),
                           "StanSensitivity/python/generate_models.py")
model_name <- GenerateSensitivityFromModel(
  base_model_name, python_script=python_script)

##################################
# Compile and run the base model.

model <- stan_model(paste(model_name, "_generated.stan", sep=""))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(model_name, "data.R", sep="."), local=stan_data)
stan_data <- as.list(stan_data)

# For now, you must use chains=1 for now to avoid confusion around get_inits.
# The script currently assumes the same number of warm-up draws as final samples.
mcmc_time <- Sys.time()
sampling_result <- sampling(
  model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples))
mcmc_time <- Sys.time() - mcmc_time
print(summary(sampling_result))

##################################
# Get the sensitivity model and sensitivity.

draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
stan_sensitivity_model <- GetStanSensitivityModel(model_name, stan_data)
sens_time <- Sys.time()
sens_result <- GetStanSensitivityFromModelFit(
  sampling_result, draws_mat, stan_sensitivity_model)
sens_time <- Sys.time()- sens_time


##################################
# Inspect the results.

# Warning: the uncertainty estimates on the sensitivity are currently
# underestimated, as they do not take into account autocorrelation in the
# MCMC chain.
tidy_results <- GetTidyResult(draws_mat, sens_result)

ggplot(filter(tidy_results$sens_norm_df, !grepl("weight", hyperparameter))) +
  geom_bar(aes(x=parameter, y=mean_sensitivity, fill=hyperparameter),
           stat="identity", position="dodge") +
  geom_errorbar(aes(x=parameter, ymin=lower_sensitivity,
                    ymax=upper_sensitivity, group=hyperparameter),
                position=position_dodge(0.9), width=0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Sensitivity / standard deviation")


ggplot(filter(tidy_results$sens_df, !grepl("weight", hyperparameter))) +
  geom_bar(aes(x=parameter, y=mean_sensitivity, fill=hyperparameter),
           stat="identity", position="dodge") +
  geom_errorbar(aes(x=parameter, ymin=lower_sensitivity,
                    ymax=upper_sensitivity, group=hyperparameter),
                position=position_dodge(0.9), width=0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
