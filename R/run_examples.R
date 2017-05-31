library(rstan)

rstan_options(auto_write = TRUE)

# Set this to be the appropriate location of the repo on your computer.
example_directory <- file.path(Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/example_models")
script_directory <- file.path(Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/R")
source(file.path(script_directory, "stan_sensitivity_lib.R"))

# Choose the basename of a model that's already been generated.  For details on how to
# generate a sensitivity model from an existing stan model, see the README in this repo.

#model_name <- file.path(example_directory, "normal_censored/normal_censored")
model_name <- file.path(example_directory, "negative_binomial/negative_binomial")

# Compile the base model.
model <- stan_model(paste(model_name, "_generated.stan", sep=""))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(model_name, "data.R", sep="."), local=stan_data)

# Get the posterior draws.
# Use chains=1 for now to avoid confusion around get_inits.  The script currently assumes
# the same number of warm-up draws as final samples.
num_samples <- 1000
result <- sampling(model, data=stan_data, chains=1, iter=num_samples * 2)
print(result)
draws_mat <- extract(result, permute=FALSE)[,1,]

stan_sensitivity_list <- GetStanSensitivityModel(model_name, stan_data)
sens_mat <- GetStanSensitivityFromModelFit(draws_mat, stan_sensitivity_list)$sens_mat

if (FALSE) {
  # Wrapping in a comment block because your analysis may not have all these variables.

  # Plot the sensitivity against whatever data value seems relevant (it's not always y).
  weight_rows <- grepl("weights", sens_param_names)
  plot(stan_data$y, sens_mat[weight_rows, 1])

  print(result)
  print(sens_mat[!weight_rows, ])
}

