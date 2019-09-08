library(devtools)
library(rstan)
rstan_options(auto_write=TRUE)

model_filename <- package_file("inst/examples/my_model.stan")
stan_data_filename <- package_file("inst/examples/my_model.data.R")

# Parse the hyperparameters block of `my_model.stan` into models that
# Stan can actuall compile.
model_name <- GenerateSensitivityFromModel(model_filename)
model <- stan_model(GetSamplingModelFilename(model_name))

# Get the original samples.
stan_data <- LoadStanData(stan_data_filename)
sampling_result <- sampling(model, data=stan_data, chains=1, iter=1000)

# Get the sensitivity result.
stan_sensitivity_list <- GetStanSensitivityModel(model_name, stan_data)
sens_result <- GetStanSensitivityFromModelFit(
    sampling_result, stan_sensitivity_list)

print(sens_result$sens_mat)
print(sens_result$sens_mat_normalized)
