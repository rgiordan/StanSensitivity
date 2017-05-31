library(rstan)

rstan_options(auto_write = TRUE)

# Set this to be the appropriate location of the repo on your computer.
example_directory <- file.path(Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/example_models")

# Choose the basename of a model that's already been generated.  For details on how to
# generate a sensitivity model from an existing stan model, see the README in this repo.

#model_name <- file.path(example_directory, "normal_censored/normal_censored")
model_name <- file.path(example_directory, "negative_binomial/negative_binomial")

# Compile the base model.
model <- stan_model(paste(model_name, "_generated.stan", sep=""))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(model_name, "data.R", sep="."), local=stan_data)

# Get the posterior draws.
# Use chains=1 for now to avoid confusion around get_inits.  The script currently assumes
# the same number of warm-up draws as final samples.
num_samples <- 1000
result <- sampling(model, data=stan_data, chains=1, iter=num_samples * 2)
print(result)

# Compile the sensitivity model and get a stanfit object.
model_sens <- stan_model(paste(model_name, "_sensitivity.stan", sep=""))
model_sens_fit <- stan(paste(model_name, "_sensitivity.stan", sep=""),
                       data=stan_data, algorithm="Fixed_param", iter=1, chains=1)


# Get the sensitivity parameters in list form.
sens_par_list <- get_inits(model_sens_fit)[[1]]
for (par in names(sens_par_list)) {
  if (par %in% names(stan_data)) {
    cat("Copying hyperparameter ", par, " from the data.\n")
    sens_par_list[[par]] <- stan_data[[par]]
  }
}

# These names help sort through the vectors of sensitivity.
param_names <- result@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)
sens_param_names <- model_sens_fit@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)

# Get the model gradients with respect to the hyperparameters (and parameters).
grad_mat <- matrix(NA, num_samples, length(sens_param_names))
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
  grad_mat[n, ] <- grad_log_prob(model_sens_fit, pars_free)
}
close(prog_bar)

# Calculate the sensitivity.
draws_mat <- extract(result, permute=FALSE)[,1,]
sens_mat <- cov(grad_mat, draws_mat)
rownames(sens_mat) <- sens_param_names

# Stan takes gradients with respect to everything in the parameters block, not just the hyperparameters.
# Remove the rows not corresponding to hyperparameters.
sens_mat <- sens_mat[setdiff(sens_param_names, param_names), ]

if (FALSE) {
  # Wrapping in a comment block because your analysis may not have all these variables.

  # Plot the sensitivity against whatever data value seems relevant (it's not always y).
  weight_rows <- grepl("weights", sens_param_names)
  plot(stan_data$y, sens_mat[weight_rows, 1])

  print(result)
  print(sens_mat[!weight_rows, ])
}

