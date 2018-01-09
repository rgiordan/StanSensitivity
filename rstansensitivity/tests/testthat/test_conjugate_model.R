library(testthat)
library(rstansensitivity)
library(rstan)
rstan_options(auto_write=TRUE)

context("rstansensitivity")

GenerateTestModels <- function() {
  base_model <- "
  data {
    real y;
  }
  hyperparameters {
    real prior_mean;
  }
  parameters {
    real mu;
  }
  model {
    mu ~ normal(prior_mean, 1.0);
    y ~ normal(mu, 1.0);
  }

  "
  
  model_name <- "/tmp/rstansensitivity_test"
  base_model_name <- paste(model_name, "stan", sep=".")
  model_file <- file(base_model_name, "w")
  cat(base_model, file=model_file)
  close(model_file)

  model_name <- GenerateSensitivityFromModel(base_model_name)

  data_text <- "
  y <- 3.0
  prior_mean <- 0.0
  "
  data_file <- file(paste(model_name, "data.R", sep="."), "w")
  cat(data_text, file=data_file)
  close(data_file)

  return(model_name)
}


test_that("conjugate_model_works", {
  test_dir <- getwd()
  model_name <- GenerateTestModels()

  model <- stan_model(paste(model_name, "_generated.stan", sep=""))

  # Load the data and hyperparameters.
  stan_data <- new.env()
  source(paste(model_name, "data.R", sep="."), local=stan_data)
  stan_data <- as.list(stan_data)

  # For now, you must use chains=1 for now to avoid confusion around get_inits.
  # The script currently assumes the same number of warm-up draws as final samples.
  iter <- 2000
  num_chains <- 3
  result <- sampling(model, data=stan_data, chains=num_chains, iter=iter)

  result_summary <- rstan::summary(result)
  mu_summary <- result_summary$summary["mu", ]

  num_samples <- (result@sim$iter - result@sim$warmup ) * num_chains
  post_var <- 1 / 2.0
  post_sd <- sqrt(post_var)
  post_se <- post_sd / sqrt(num_samples)
  post_mean <- 0.5 * (stan_data$prior_mean + mean(stan_data$y))

  # Sanity checks.
  expect_true(abs(post_mean - mu_summary["mean"]) / (3 * post_se)  < 1)
  expect_true(abs(post_sd - mu_summary["sd"]) / (3 * post_se)  < 1)
  
  # Check the sensitivity.
  stan_sensitivity_list <- GetStanSensitivityModel(model_name, stan_data)
  sens_result <- GetStanSensitivityFromModelFit(result, stan_sensitivity_list)
  sens_mat <- sens_result$sens_mat
  sens_mat_normalized <- sens_result$sens_mat_normalized

  mean_sens <- sens_mat["prior_mean", "mu"]
  expect_true(abs(post_var - mean_sens) / (3 * post_se) < 1)
})
