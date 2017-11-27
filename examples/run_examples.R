library(rstan)
library(rstansensitivity)
library(ggplot2)
library(gridExtra)

rstan_options(auto_write=TRUE)

# Set this to be the appropriate location of the repo on your computer.
example_directory <- file.path(
  Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/examples/example_models")

model_name <- file.path(example_directory, "negative_binomial/negative_binomial")
num_warmup_samples <- 5000
num_samples <- 5000

##################################
# Compile and run the base model.

model <- stan_model(paste(model_name, "_generated.stan", sep=""))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(model_name, "data.R", sep="."), local=stan_data)
stan_data <- as.list(stan_data)

# For now, you must use chains=1 for now to avoid confusion with stan's get_inits.
mcmc_time <- Sys.time()
sampling_result <- sampling(
  model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples))
mcmc_time <- Sys.time() - mcmc_time
print(summary(sampling_result))


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

