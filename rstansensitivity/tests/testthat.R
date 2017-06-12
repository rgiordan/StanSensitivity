library(testthat)
library(rstansensitivity)
library(rstan)
rstan_options(auto_write=TRUE)

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

base_model_generated <- "
data {
  real y;

  // Hyperparameters:
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

base_model_sensitivity <- "
data {
  real y;
}
parameters {
  real mu;

  // Hyperparameters:
  real prior_mean;
}
model {
  mu ~ normal(prior_mean, 1.0);
  y ~ normal(mu, 1.0);
}
"

# To generate the files:
# python_script <- file.path(Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/python/generate_models.py")
# model_name <- GenerateSensitivityFromModel(base_model_name, python_script=python_script)

model_name <- "/tmp/rstansensitivity_test"
base_model_name <- paste(model_name, "stan", sep=".")
model_file <- file(base_model_name, "w")
cat(base_model, file=model_file)
close(model_file)

model_file <- file(paste(model_name, "generated.stan", sep="_"), "w")
cat(base_model_generated, file=model_file)
close(model_file)

model_file <- file(paste(model_name, "sensitivity.stan", sep="_"), "w")
cat(base_model_sensitivity, file=model_file)
close(model_file)


data_text <- "
y <- 3.0
prior_mean <- 0.0
"
data_file <- file(paste(model_name, "data.R", sep="."), "w")
cat(data_text, file=data_file)
close(data_file)

#test_check("interpret_model", {

model <- stan_model(paste(model_name, "_generated.stan", sep=""))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(model_name, "data.R", sep="."), local=stan_data)
stan_data <- as.list(stan_data)

# For now, you must use chains=1 for now to avoid confusion around get_inits.
# The script currently assumes the same number of warm-up draws as final samples.
num_warmup_samples <- 1000
num_samples <- 1000
result <- sampling(model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples))

draws_mat <- extract(result, permute=FALSE)[,1,]

post_var <- 1 / 2.0
post_sd <- sqrt(post_var)
post_se <- post_sd / sqrt(num_samples)
post_mean <- 0.5 * (stan_data$prior_mean + mean(stan_data$y))

# Sanity checks.
expect_equal(post_mean, mean(draws_mat[, "mu"]), tolerance=3 * post_se)
expect_equal(post_sd, sd(draws_mat[, "mu"]), tolerance=3 * post_se)

# Check the sensitivity.
stan_sensitivity_list <- GetStanSensitivityModel(model_name, stan_data)
sens_result <- GetStanSensitivityFromModelFit(
  draws_mat, stan_sensitivity_list, num_warmup_samples=num_warmup_samples)
sens_mat <- sens_result$sens_mat
sens_mat_normalized <- sens_result$sens_mat_normalized

sens_mat[, "mu", drop=FALSE]

#})




model_sens_fit <- stan_sensitivity_list$model_sens_fit
param_names <- stan_sensitivity_list$param_names
sens_param_names <- stan_sensitivity_list$sens_param_names
sens_par_list <- stan_sensitivity_list$sens_par_list

# Get the model gradients with respect to the hyperparameters (and parameters).
num_samples <- nrow(draws_mat)
grad_mat <- matrix(NA, num_samples, length(sens_param_names))
lp_vec <- rep(NA, num_samples)
cat("Evaluating log gradients at the MCMC draws.\n")
prog_bar <- txtProgressBar(min=1, max=num_samples, style=3)
for (n in 1:num_samples) {
  setTxtProgressBar(prog_bar, value=n)
  par_list <- get_inits(result, iter=n + num_samples)[[1]]
  for (par in ls(par_list)) {
    # Note that get_inits is currently broken:
    # https://github.com/stan-dev/rstan/issues/417
    # ...but this has seemed to fix it (so far):
    sens_par_list[[par]] <- as.numeric(par_list[[par]])
  }
  pars_free <- unconstrain_pars(model_sens_fit, sens_par_list)
  glp <- grad_log_prob(model_sens_fit, pars_free)
  grad_mat[n, ] <- glp
  lp_vec[n] <- attr(glp, "log_prob")
}
close(prog_bar)

sens_mat <- cov(grad_mat, draws_mat)
rownames(sens_mat) <- sens_param_names
sens_mat <- sens_mat[setdiff(sens_param_names, param_names), , drop=FALSE]


