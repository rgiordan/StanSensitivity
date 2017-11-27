# This script demonstrates the basic usage of the rstansensitivity library
# using the negative_binomial example.

library(rstan)
library(rstansensitivity)
library(ggplot2)
library(gridExtra)

rstan_options(auto_write=TRUE)

# Run from anywhere in the StanSensitivity repository.
git_repo <- system("git rev-parse --show-toplevel", intern=TRUE)
example_dir <- file.path(git_repo, "examples/normal_censored/")
model_name <- file.path(example_dir, "models/normal_censored")
num_warmup_samples <- 10000
num_samples <- 10000

##################################
# Compile and run the base model.

model <- stan_model(paste(model_name, "_generated.stan", sep=""))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(example_dir, "normal_censored.data.R", sep=""), local=stan_data)
stan_data <- as.list(stan_data)

# For now, rstansensitivity only supports one chain.
sampling_result <- sampling(
  model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples))
print(sampling_result)

stan_sensitivity_model <- GetStanSensitivityModel(model_name, stan_data)
sens_result <- GetStanSensitivityFromModelFit(sampling_result, stan_sensitivity_model)
tidy_results <- GetTidyResult(sens_result)

# There's considerable sensitivity with repsect to y_var.
filter(tidy_results, !grepl("weights", hyperparameter))
stan_data$y_var

# Look at the pattern of data weight sensitivities.
weight_sens <-
  filter(tidy_results, grepl("weights", hyperparameter)) %>%
  mutate(row=as.integer(gsub("weights\\.", "", hyperparameter))) %>%
  inner_join(data.frame(row=1:length(stan_data$y), y=stan_data$y), by="row")

ggplot(weight_sens) +
  geom_point(aes(x=y, y=normalized_sensitivity)) +
  geom_errorbar(aes(x=y,
                    ymin=normalized_sensitivity - 2 * normalized_sensitivity_se,
                    ymax=normalized_sensitivity + 2 * normalized_sensitivity_se))
