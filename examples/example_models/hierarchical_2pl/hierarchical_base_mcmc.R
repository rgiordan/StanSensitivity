library(rstan)
rstan_options(auto_write=TRUE)
library(ggplot2)
library(dplyr)
library(reshape2)
library(rstansensitivity)

library(jsonlite)


# Set this to be the appropriate location of the repo on your computer.
example_directory <- file.path(
  Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/examples/example_models/hierarchical_2pl")

mcmc_results_json_filename <- file.path(example_directory, "hierarchical_2pl_mcmc_results.json")

library(mirt)
sat <- key2binary(SAT12,
                  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))

stan_data <- list(I = ncol(sat),
                  J = nrow(sat),
                  N = length(sat),
                  ii = rep(1:ncol(sat), each = nrow(sat)),
                  jj = rep(1:nrow(sat), times = ncol(sat)),
                  y = as.vector(sat))

# Set the hyperparameters.
stan_data$lkj_concentration <- 4.0
stan_data$mu_loc <- 0.
stan_data$mu_1_scale <-1.
stan_data$mu_2_scale <- 5.
stan_data$tau_loc <- 0.1
stan_data$theta_loc <- 0.
stan_data$theta_scale <- 1.

base_model_name <- file.path(example_directory, "hierarchical_2pl.stan")
python_script <- file.path(Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/python/generate_models.py")
model_name <- GenerateSensitivityFromModel(base_model_name, python_script=python_script)
model <- stan_model(paste(model_name, "_generated.stan", sep=""))

num_warmup_samples <- num_samples <- 1000
mcmc_time <- Sys.time()
sampling_result <- sampling(model, data=stan_data, chains=1, iter=num_warmup_samples + num_samples)
mcmc_time <- Sys.time() - mcmc_time

# Write to JSON so python can read it.

sampling_result_summary <- data.frame(summary(sampling_result)$summary)
sampling_result_json <- list(pars=rownames(sampling_result_summary),
                             mean=sampling_result_summary$mean,
                             sd=sampling_result_summary$sd,
                             se_mean=sampling_result_summary$se_mean)

json_file <- file(mcmc_results_json_filename, "w")
json_list <- toJSON(list(mcmc_results=sampling_result_json, stan_data=stan_data))
write(json_list, file=json_file)
close(json_file)

