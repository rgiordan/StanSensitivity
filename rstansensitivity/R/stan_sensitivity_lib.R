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
GetStanSensitivityModel <- function(sampling_result, model_name, stan_data) {
    # Use a "model" with no model block to get a valid parameter list with
    # get_inits.  This way there is no worry about invalid intializations.
    model_sens_params <-
        stan(paste(model_name, "_sensitivity_parameters.stan", sep=""),
                   data=stan_data, algorithm="Fixed_param",
                   iter=1, chains=1)

    sens_par_list <- get_inits(model_sens_params)[[1]]
    model_par_list <- get_inits(sampling_result)[[1]]

    # Get the sensitivity parameters in list form.
    for (par in names(model_par_list)) {
      cat("Copying parameter '", par, "' from the sampler\n", sep="")
      sens_par_list[[par]] <- model_par_list[[par]]
    }
    for (par in names(sens_par_list)) {
        if (par %in% names(stan_data)) {
            cat("Copying hyperparameter '", par,
                 "' from the data block.\n", sep="")
            sens_par_list[[par]] <- stan_data[[par]]
        }
    }

    model_sens_fit <- stan(paste(model_name, "_sensitivity.stan", sep=""),
                           data=stan_data, algorithm="Fixed_param",
                           iter=1, chains=1, init=list(sens_par_list))

    # These names help sort through the vectors of sensitivity.
    param_names <-
        sampling_result@.MISC$stan_fit_instance$unconstrained_param_names(
            FALSE, FALSE)
    sens_param_names <-
        model_sens_fit@.MISC$stan_fit_instance$unconstrained_param_names(
            FALSE, FALSE)

    return(list(model_sens_fit=model_sens_fit,
                model_sens_params=model_sens_params,
                param_names=param_names,
                sens_param_names=sens_param_names,
                sens_par_list=sens_par_list))
}


# Evaluate at draws using the hyperparameters in sens_par_list.
EvaluateAtDraws <- function(
    sampling_result, draws_mat, stan_sensitivity_list, sens_par_list,
    compute_grads=FALSE) {

  num_warmup_samples <- attr(sampling_result, "sim")$warmup
  model_sens_fit <- stan_sensitivity_list$model_sens_fit
  sens_param_names <- stan_sensitivity_list$sens_param_names

  # Get the model gradients with respect to the hyperparameters (and parameters).
  num_samples <- nrow(draws_mat)
  lp_vec <- rep(NA, num_samples)
  if (compute_grads) {
    grad_mat <- matrix(NA, length(sens_param_names), num_samples)
  } else {
    grad_mat <- matrix()
  }

  cat("Evaluating model at the MCMC draws.\n")
  prog_bar <- txtProgressBar(min=1, max=num_samples, style=3)
  for (n in 1:num_samples) {
    setTxtProgressBar(prog_bar, value=n)
    par_list <- get_inits(sampling_result, iter=n + num_warmup_samples)[[1]]
    for (par in ls(par_list)) {
      # Note that get_inits is currently broken:
      # https://github.com/stan-dev/rstan/issues/417
      # ...but this has seemed to fix it (so far):
      if (length(dim(sens_par_list[[par]])) >= 2) {
        sens_par_list[[par]] <-
          array(unlist(par_list[[par]]), dim(par_list[[par]]))
      } else {
        sens_par_list[[par]] <- as.numeric(par_list[[par]])
      }
    }
    pars_free <- unconstrain_pars(model_sens_fit, sens_par_list)

    if (compute_grads) {
      glp <- grad_log_prob(model_sens_fit, pars_free)
      grad_mat[, n] <- glp
      lp_vec[n] <- attr(glp, "log_prob")
    } else {
      lp_vec[n] <- log_prob(model_sens_fit, pars_free)
    }
  }
  close(prog_bar)
  return(list(lp_vec=lp_vec, grad_mat=grad_mat))
}


# Process the results of GetStanSensitivityModel into a sensitivity matrix.
GetStanSensitivityFromModelFit <- function(
    sampling_result, draws_mat, stan_sensitivity_list) {

    # Get the model gradients with respect to the hyperparameters (and parameters).
    model_at_draws <- EvaluateAtDraws(
            sampling_result, draws_mat, stan_sensitivity_list,
            stan_sensitivity_list$sens_par_list, compute_grads=TRUE)

    param_names <- stan_sensitivity_list$param_names
    sens_param_names <- stan_sensitivity_list$sens_param_names

    # Stan takes gradients with respect to everything in the parameters block,
    # not just the hyperparameters.  Remove the rows not corresponding to
    # hyperparameters.
    grad_mat <- model_at_draws$grad_mat
    rownames(grad_mat) <- sens_param_names
    grad_mat <- grad_mat[setdiff(sens_param_names, param_names),, drop=FALSE]

    # Calculate the sensitivity.
    sens_mat <- GetSensitivityFromGrads(grad_mat, draws_mat)

    # Normalize by the marginal standard deviation.
    sens_mat_normalized <- NormalizeSensitivityMatrix(sens_mat, draws_mat)

    return(list(sens_mat=sens_mat, sens_mat_normalized=sens_mat_normalized,
                grad_mat=grad_mat, lp_vec=model_at_draws$lp_vec))
}


GetImportanceSamplingFromModelFit <- function(
    sampling_result, draws_mat, stan_sensitivity_list,
    imp_sens_par_list, lp_vec=NULL) {

    if (is.null(lp_vec)) {
        lp_vec <- EvaluateAtDraws(
            sampling_result, draws_mat, stan_sensitivity_list,
            stan_sensitivity_list$sens_par_list, compute_grads=FALSE)$lp_vec
    }
    # Get the model gradients with respect to the hyperparameters (and parameters).
    imp_lp_vec <- EvaluateAtDraws(
            sampling_result, draws_mat, stan_sensitivity_list,
            imp_sens_par_list, compute_grads=FALSE)$lp_vec

    imp_weights <- exp(imp_lp_vec - lp_vec)
    imp_weights <- imp_weights / sum(imp_weights)

    eff_num_imp_samples <- 1 / sum(imp_weights ^ 2)

    imp_means <- colSums(imp_weights * draws_mat)

    return(list(imp_weights=imp_weights, imp_lp_vec=imp_lp_vec, lp_vec=lp_vec,
                eff_num_imp_samples=eff_num_imp_samples,
                imp_means=imp_means))
}


#########################################
# Bootstrap the covariance calculations.

GetSensitivityFromGrads <- function(grad_mat, draws_mat) {
  # This should in fact match cov() but without having to transpose,
  # which gives a speedup.

  # Should match:
  # sens_mat <- cov(t(grad_mat), draws_mat)

  grad_means <- rowMeans(grad_mat)
  draw_means <- colMeans(draws_mat)

  n <- nrow(draws_mat)
  sens_mat <- (grad_mat %*% draws_mat) / (n - 1) -
               grad_means %*% t(draw_means) * n / (n - 1)

  rownames(sens_mat) <- rownames(grad_mat)
  colnames(sens_mat) <- colnames(draws_mat)

  return(sens_mat)
}


NormalizeSensitivityMatrix <- function(sens_mat, draws_mat) {
  draw_sds <- sqrt(apply(draws_mat, sd, MARGIN=2))
  return(sens_mat / rep(draw_sds, each=nrow(sens_mat)))
}


BootstrapSensitivityMatrix <- function(draws_mat, grad_mat, alpha=0.1, num_boot=200) {
  cat("Bootstrapping sensitivity matrix.")
  prog_bar <- txtProgressBar(min=1, max=num_boot, style=3)
  num_boot <- 200
  sens_mat_dim <-
  sens_mat_array <- array(NA, dim=c(num_boot, dim(sens_mat)))
  sens_mat_normalized_array <- array(NA, dim=c(num_boot, dim(sens_mat)))
  for (boot in 1:num_boot) {
    setTxtProgressBar(prog_bar, value=boot)
    w <- sample.int(n=nrow(draws_mat), size=nrow(draws_mat), replace=TRUE)
    sens_mat_boot <-
      GetSensitivityFromGrads(grad_mat[, w, drop=FALSE],
                              draws_mat[w, , drop=FALSE])
    sens_mat_array[boot,,] <- sens_mat_boot
    sens_mat_normalized_array[boot,,] <-
      NormalizeSensitivityMatrix(sens_mat_boot, draws_mat[w, , drop=FALSE])
  }
  close(prog_bar)
  return(list(sens_mat_array=sens_mat_array,
              sens_mat_normalized_array=sens_mat_normalized_array))
}


GetArrayQuantiles <- function(sens_mat_array, alpha=0.1) {
  lower <- apply(sens_mat_array, MARGIN=c(2, 3),
                 function(x) { quantile(x, alpha) })
  upper <- apply(sens_mat_array, MARGIN=c(2, 3),
                 function(x) { quantile(x, 1 - alpha) })
  return(list(lower=lower, upper=upper))
}
