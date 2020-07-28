library(testthat)
library(rstansensitivity)
library(rstan)
rstan_options(auto_write=TRUE)

context("rstansensitivity")

CheckTidyResults <- function(tidy_results) {
    expect_equal(tidy_results$sensitivity_upper,
        tidy_results$sensitivity + 2 * tidy_results$sensitivity_se,
        tolerance=1e-8)
    expect_equal(tidy_results$sensitivity_lower,
        tidy_results$sensitivity - 2 * tidy_results$sensitivity_se,
        tolerance=1e-8)
    expect_equal(tidy_results$normalized_sensitivity_upper,
        tidy_results$normalized_sensitivity +
        2 * tidy_results$normalized_sensitivity_se,
        tolerance=1e-8)
    expect_equal(tidy_results$normalized_sensitivity_lower,
        tidy_results$normalized_sensitivity -
        2 * tidy_results$normalized_sensitivity_se,
        tolerance=1e-8)
    PlotSensitivities(tidy_results)
}


################################################
################################################
################################################
test_that("basic_model_works", {

# A silly model just for testing.
base_model <- "
data {
int N;
int K;
real y[N];
real x[K];
}
hyperparameters {
real prior_mean[K];
real prior_var[K];
}
parameters {
real beta[5];
}
model {
beta ~ normal(prior_mean, prior_var);
y ~ normal(dot_product(beta, x), 1.0);
}

"

model_name <- "/tmp/rstansensitivity_post_test"
base_model_name <- paste(model_name, "stan", sep=".")
model_file <- file(base_model_name, "w")
cat(base_model, file=model_file)
close(model_file)

model_name <- GenerateSensitivityFromModel(base_model_name)
model <- stan_model(GetSamplingModelFilename(model_name))

set.seed(42)
N <- 10
K <- 5
stan_data <- list(
    N=N, K=K, y=runif(N), x=runif(K),
    prior_mean=runif(K), prior_var=exp(runif(K))
)

iter <- 500
num_chains <- 2
sampling_result <- sampling(model, data=stan_data, chains=num_chains, iter=iter)

# Mostly I will just check that these run succesfully.
stan_sensitivity_list <- GetStanSensitivityModel(model_name, stan_data)
sens_result <- GetStanSensitivityFromModelFit(sampling_result, stan_sensitivity_list)
tidy_results <- GetTidyResult(sens_result)
CheckTidyResults(tidy_results)
hyperparam_df <- GetHyperparameterDataFrame(stan_sensitivity_list, stan_data)
expect_equal(
    filter(hyperparam_df,
        grepl("^prior_mean", hyperparameter))$hyperparameter_val,
    stan_data$prior_mean)
expect_equal(
    filter(hyperparam_df,
        grepl("^prior_var", hyperparameter))$hyperparameter_val,
    stan_data$prior_var)

expect_equivalent(as.character(unique(hyperparam_df$hyperparameter)),
                  as.character(unique(tidy_results$hyperparameter)))


# Check the transform.
transform_mat <- matrix(0, 2, 2 * K)
colnames(transform_mat) <- rownames(sens_result$sens_mat)
transform_mat[, "prior_mean.1"] <- 1
transform_mat[2, "prior_mean.2"] <- 1
rownames(transform_mat) <- c("prior_mean.1.only", "prior_mean.1p2")
trans_sens_result <- TransformSensitivityResult(sens_result, transform_mat)
trans_tidy_results <- GetTidyResult(trans_sens_result)
numeric_cols <- paste("sensitivity", c("", "_se", "_upper", "_lower"), sep="")
numeric_cols <- c(numeric_cols, paste("normalized", numeric_cols, sep="_"))
expect_equivalent(
    filter(trans_tidy_results,
        hyperparameter == "prior_mean.1.only")[, numeric_cols],
    filter(tidy_results,
        hyperparameter == "prior_mean.1")[, numeric_cols])
sens_cols <- c("sensitivity", "normalized_sensitivity")
expect_equivalent(
    filter(trans_tidy_results,
        hyperparameter == "prior_mean.1p2")[, sens_cols],
    filter(tidy_results, hyperparameter == "prior_mean.1")[, sens_cols] +
    filter(tidy_results, hyperparameter == "prior_mean.2")[, sens_cols])

# Check the MCMC dataframe.
mcmc_df <- GetMCMCDataFrame(sampling_result)
expect_equivalent(as.character(unique(mcmc_df$parameter)),
                  as.character(unique(tidy_results$parameter)))

# Check the perturbation.
stan_data_perturb <- stan_data
stan_data_perturb$prior_mean[1] <- stan_data_perturb$prior_mean[1] + 1
sens_pred_df <- PredictSensitivityFromStanData(
    stan_sensitivity_list, sens_result,
    stan_data, stan_data_perturb, "prior_mean.1.only")
expect_equivalent(
    filter(sens_pred_df,
        hyperparameter == "prior_mean.1.only")[, numeric_cols],
    filter(tidy_results,
        hyperparameter == "prior_mean.1")[, numeric_cols])

})
