# This script demonstrates the basic usage of the rstansensitivity library
# using the negative_binomial example.

library(rstan)
library(rstansensitivity)
library(ggplot2)
library(gridExtra)

rstan_options(auto_write=TRUE)
options(warn=1) # Display Stan warnings as they occur.

# Run from anywhere in the StanSensitivity repository.
git_repo <- system("git rev-parse --show-toplevel", intern=TRUE)
example_dir <- file.path(git_repo, "examples/normal_censored/")
model_name <- file.path(example_dir, "models/normal_censored")
iters <- 20000

##################################
# Compile and run the base model.

model <- stan_model(paste(model_name, "_generated.stan", sep=""))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(example_dir, "normal_censored.data.R", sep=""), local=stan_data)
stan_data <- as.list(stan_data)

# For now, rstansensitivity only supports one chain.
sampling_result <- sampling(
  model, data=stan_data, chains=1, iter=iters)
print(sampling_result)

stan_sensitivity_model <- GetStanSensitivityModel(model_name, stan_data)
sens_result <- GetStanSensitivityFromModelFit(sampling_result, stan_sensitivity_model)
tidy_results <- GetTidyResult(sens_result)
PlotSensitivities(filter(tidy_results, !grepl("weights", hyperparameter)))

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


###############################
# Refit to check the linearity assumption

StanResultToDataframe <- function(sampling_result) {
  sampling_df <- data.frame(summary(sampling_result)$summary)
  sampling_df$parameter <- make.names(rownames(sampling_df))
  return(sampling_df)
}

hyperparam_name <- "y_var"
hyperparam_min_offset <- -0.8
hyperparam_max_offset <- 1.0
num_refits <- 20

hyperparam_vals <- seq(stan_data[[hyperparam_name]] + hyperparam_min_offset,
                       stan_data[[hyperparam_name]] + hyperparam_max_offset,
                       length.out=num_refits)

sens_params <- unique(tidy_results$parameter)

summary_list <- list()
sampling_result_list <- list()
refit_time <- Sys.time()
for (val in hyperparam_vals) {
  cat(sprintf("Value: %f", val))
  stan_data_perturbed <- stan_data
  stan_data_perturbed[[hyperparam_name]] <- val
  sampling_result_perturbed <-
    sampling(model, data=stan_data_perturbed, iter=iters, chains=1)
  sampling_result_list[[length(sampling_result_list) + 1]] <- sampling_result_perturbed
  sens_summary_peturbed <-
    StanResultToDataframe(sampling_result_perturbed) %>%
    mutate(hyperparameter=hyperparam_name, hyperparameter_val=val)
  summary_list[[length(summary_list) + 1]] <- sens_summary_peturbed
}
refit_time <- Sys.time() - refit_time

perturbed_df <-
  do.call(rbind, summary_list) %>%
  filter(parameter %in% sens_params) %>%
  inner_join(tidy_results, by=c("parameter", "hyperparameter")) %>%
  inner_join(StanResultToDataframe(sampling_result),
             by="parameter", suffix=c("", ".orig"))

print(ggplot(perturbed_df) +
        geom_point(aes(x=hyperparameter_val, y=mean - mean.orig, color=parameter, group=parameter)) +
        geom_errorbar(aes(x=hyperparameter_val,
                          ymin=mean - mean.orig - 2 * se_mean,
                          ymax=mean - mean.orig + 2 * se_mean,
                          color=parameter, group=parameter), width=0.1) +
        geom_line(aes(x=hyperparameter_val,
                      y=(hyperparameter_val - stan_data[[hyperparam_name]]) * sensitivity,
                      color=parameter, group=parameter)) +
        geom_ribbon(aes(x=hyperparameter_val,
                        ymin=(hyperparameter_val - stan_data[[hyperparam_name]]) *
                          (sensitivity - 2 * sensitivity_se) - 2 * se_mean.orig,
                        ymax=(hyperparameter_val - stan_data[[hyperparam_name]]) *
                          (sensitivity + 2 * sensitivity_se) + 2 * se_mean.orig,
                        fill=parameter, group=parameter), color=NA, alpha=0.2) +
        facet_grid(~ parameter))

# Save the results
results_file <- sprintf("%s_results.Rdata", model_name)
save(stan_data, hyperparam_vals, hyperparam_name,
     sampling_result, tidy_results, perturbed_df, file=results_file)
