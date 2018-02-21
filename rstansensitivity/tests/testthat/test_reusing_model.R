library(testthat)
library(rstansensitivity)
library(rstan)
rstan_options(auto_write=TRUE)

context("rstansensitivity")

GenerateTestModels <- function() {
    orig_model <- "
    data {
        real y;
    }
    parameters {
        real mu;
    }
    model {
        mu ~ normal(0.0, 1.0);
        y ~ normal(mu, 1.0);
    }

    "

    sens_model <- "
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

    bad_model <- "
    data {
        real y;
        real prior_var;
    }
    parameters {
        real mu;
    }
    model {
        mu ~ normal(0.0, prior_var);
        y ~ normal(mu, 1.0);
    }

    "

    WriteModelFile <- function(model_path, model_code) {
        model_file <- file(model_path, "w")
        cat(model_code, file=model_file)
        close(model_file)        
    }

    WriteModelFile(
        orig_model_file <- "/tmp/rstansensitivity_reuse_test_orig.stan",
        orig_model)
    WriteModelFile(
        sens_model_file <- "/tmp/rstansensitivity_reuse_test_sens.stan",
        sens_model)
    WriteModelFile(
        bad_model_file <- "/tmp/rstansensitivity_reuse_test_bad.stan",
        bad_model)

    model_name <- GenerateSensitivityFromModel(sens_model_file)

    return(list(model_name=model_name,
                orig_model_file=orig_model_file,
                bad_model_file=bad_model_file))
}


test_that("reusing_model_works", {
  set.seed(42)
  test_dir <- getwd()
  model_files <- GenerateTestModels()

  # Compare the generated model to the original model in which the
  # hyperparameter was hard-encoded, and a "bad" model which differs from
  # either of the two.
  model_gen <- stan_model(GetSamplingModelFilename(model_files$model_name))
  model_orig <- stan_model(model_files$orig_model_file)
  model_bad <- stan_model(model_files$orig_model_file)

  # Load the data and hyperparameters.
  stan_data_gen <- list(y=3.0, prior_mean=0.0)
  stan_data_orig <- list(y=3.0)
  stan_data_bad <- list(y=3.0, prior_var=2.0)
  stan_data_notbad <- list(y=3.0, prior_var=1.0)

  iter <- 50
  num_chains <- 3
  sampling_result <- sampling(model_orig, data=stan_data_orig,
                     chains=num_chains, iter=iter)
        
  orig_ok <- CheckModelEquivalence(
      sampling_result, model_orig, model_gen, stan_data_orig, stan_data_gen)             
  bad_ok <- CheckModelEquivalence(
      sampling_result, model_bad, model_gen, stan_data_bad, stan_data_gen)             
  notbad_ok <- CheckModelEquivalence(
      sampling_result, model_bad, model_gen, stan_data_notbad, stan_data_gen)             
  expect(orig_ok)
  expect(!bad_ok)
  expect(notbad_ok)

  
})
