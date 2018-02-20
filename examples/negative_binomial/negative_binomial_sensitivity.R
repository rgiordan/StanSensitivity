# This script demonstrates the basic usage of the rstansensitivity library
# using the negative_binomial example.

library(rstan)
library(rstansensitivity)
library(ggplot2)
library(gridExtra)

rstan_options(auto_write=TRUE)

# Run from anywhere in the StanSensitivity repository.
git_repo <- system("git rev-parse --show-toplevel", intern=TRUE)
example_dir <- file.path(git_repo, "examples/negative_binomial/")
model_name <- GenerateSensitivityFromModel(
  file.path(example_dir, "models/negative_binomial.stan"))

##################################
# Compile and run the base model.

model <- stan_model(GetSamplingModelFilename(model_name))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(example_dir, "negative_binomial.data.R", sep=""), local=stan_data)
stan_data <- as.list(stan_data)

sampling_result <- sampling(
  model, data=stan_data, chains=3, iter=3000)
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
