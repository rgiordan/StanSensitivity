library(rstan)
library(rstansensitivity)

rstan_options(auto_write=TRUE)

# Set this to be the appropriate location of the repo on your computer.
example_directory <- file.path(
  Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/examples/example_models")

#base_model_name <- file.path(example_directory, "negative_binomial/negative_binomial.stan")
base_model_name <- file.path(example_directory, "normal_censored/normal_censored.stan")

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
num_warmup_samples <- 5000
num_samples <- 5000
sampling_result <- sampling(model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples))
print(summary(sampling_result))

##################################
# Get the sensitivity model and sensitivity.

draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
stan_sensitivity_list <- GetStanSensitivityModel(sampling_result, model_name, stan_data)
sens_result <- GetStanSensitivityFromModelFit(
  sampling_result, draws_mat, stan_sensitivity_list, num_warmup_samples=num_warmup_samples)


##################################
# Inspect the results.

sens_mat <- sens_result$sens_mat
sens_mat_normalized <- sens_result$sens_mat_normalized

# Plot the sensitivity against whatever data value seems relevant (it's not always y).
weight_rows <- grepl("weights", rownames(sens_mat))
plot(stan_data$y, sens_mat_normalized[weight_rows, 1])

# Look at the sensitivity to other hyperparameters.
print(sampling_result)
print(sens_mat[!weight_rows, , drop=FALSE])
print(sens_mat_normalized[!weight_rows, , drop=FALSE])

