library(devtools)
library(rstan)
rstan_options(auto_write=TRUE)

model_filename <- package_file("inst/examples/my_model.stan")
stan_data_filename <- package_file("inst/examples/my_model.data.R")

model_name <- GenerateSensitivityFromModel(model_filename)
model <- stan_model(GetSamplingModelFilename(model_name))

stan_data <- LoadStanData(stan_data_filename)
sampling_result <- sampling(model, data=stan_data, chains=1, draws=1000)
