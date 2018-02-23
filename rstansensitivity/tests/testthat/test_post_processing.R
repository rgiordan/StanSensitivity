library(testthat)
library(rstansensitivity)
library(rstan)
rstan_options(auto_write=TRUE)

context("rstansensitivity")

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

iter <- 100
num_chains <- 2
result <- sampling(model, data=stan_data, chains=num_chains, iter=iter)


})
