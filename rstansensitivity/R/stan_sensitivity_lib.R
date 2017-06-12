library(rstan)


GenerateSensitivityFromModel <- function(
        base_model_name,
        python_script="StanSensitivity/python/generate_models.py") {

    model_suffix <-
        substr(base_model_name, nchar(base_model_name) - 4, nchar(base_model_name))
    stopifnot(model_suffix == ".stan")
    system(paste(python_script, " --base_model=", base_model_name, sep=""))
    model_name <- sub("\\.stan$", "", base_model_name)
    return(model_name)
}


# Compile the sensitivity model and get a stanfit object and related information.
GetStanSensitivityModel <- function(model_name, stan_data) {
  model_sens <- stan_model(paste(model_name, "_sensitivity.stan", sep=""))
  model_sens_fit <- stan(paste(model_name, "_sensitivity.stan", sep=""),
                         data=stan_data, algorithm="Fixed_param",
                         iter=1, chains=1)

  # Get the sensitivity parameters in list form.
  sens_par_list <- get_inits(model_sens_fit)[[1]]
  for (par in names(sens_par_list)) {
    if (par %in% names(stan_data)) {
      cat("Copying hyperparameter '", par, "' from the data.\n", sep="")
      sens_par_list[[par]] <- stan_data[[par]]
    }
  }

  # These names help sort through the vectors of sensitivity.
  param_names <-
    result@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)
  sens_param_names <-
    model_sens_fit@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)

  return(list(model_sens_fit=model_sens_fit,
              param_names=param_names,
              sens_param_names=sens_param_names,
              sens_par_list=sens_par_list))
}


# Process the results of GetStanSensitivityModel into a sensitivity matrix.
GetStanSensitivityFromModelFit <- function(
    draws_mat, stan_sensitivity_list, num_warmup_samples=nrow(draws_mat)) {

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
    par_list <- get_inits(result, iter=n + num_warmup_samples)[[1]]
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

  # Calculate the sensitivity.
  sens_mat <- cov(grad_mat, draws_mat)
  rownames(sens_mat) <- sens_param_names

  # Stan takes gradients with respect to everything in the parameters block,
  # not just the hyperparameters.  Remove the rows not corresponding to
  # hyperparameters.
  sens_mat <- sens_mat[setdiff(sens_param_names, param_names), ]

  # Normalize by the marginal standard deviation.
  draws_sd <- sapply(1:ncol(draws_mat), function(col) sd(draws_mat[, col]))
  sens_mat_normalized <- t(t(sens_mat) / draws_sd)

  return(list(sens_mat=sens_mat, sens_mat_normalized=sens_mat_normalized,
              grad_mat=grad_mat, lp_vec=lp_vec))
}
