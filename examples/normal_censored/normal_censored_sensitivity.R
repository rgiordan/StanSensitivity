# This script demonstrates the basic usage of the rstansensitivity library
# using the negative_binomial example.

library(rstan)
library(rstansensitivity)
library(ggplot2)
library(gridExtra)

rstan_options(auto_write=TRUE)

# Run from anywhere in the StanSensitivity repository.
git_repo <- system("git rev-parse --show-toplevel", intern=TRUE)
example_directory <- file.path(git_repo, "examples/example_models")

model_name <- file.path(example_directory, "negative_binomial/negative_binomial")
num_warmup_samples <- 10000
num_samples <- 10000

##################################
# Compile and run the base model.

model <- stan_model(paste(model_name, "_generated.stan", sep=""))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(model_name, "data.R", sep="."), local=stan_data)
stan_data <- as.list(stan_data)

# For now, rstansensitivity only supports one chain.
sampling_result <- sampling(
  model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples))
print(sampling_result)


##################################
# Get the sensitivity model and sensitivity.

stan_sensitivity_model <- GetStanSensitivityModel(model_name, stan_data)
sens_result <- GetStanSensitivityFromModelFit(sampling_result, stan_sensitivity_model)
tidy_results <- GetTidyResult(sens_result)

grid.arrange(
  PlotSensitivities(tidy_results, normalized=TRUE) +
    ggtitle("Normalized sensitivities")
  ,
  PlotSensitivities(tidy_results, normalized=FALSE) +
    ggtitle("Sensitivities")
  ,
  ncol=2)
