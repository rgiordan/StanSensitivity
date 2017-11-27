library(rstan)
library(dplyr)
library(reshape2)


#' Generate stan models for sensitivity calculations from a model with a
#' hyperparameters block.
#'
#' @param base_model_name The name of the model with a hyperparameters block,
#' including the .stan suffix.
#' @param python_script The location of the generate_models.py from this
#' repository.
#' @return The \code{model_name} which can be passed as an argument into other
#' functions in this library.  The function also generates the sensitivity model
#' files in the same location as the original base model.
#' @export
#' @examples
#' GenerateSensitivityFromModel(
#'     "models/my_model.stan",
#'     "~/git_repos/StanSensitivity/python/generate_models.py")
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


#' Compile the sensitivity model and get a stanfit object and related
#' information.
#'
#' @param model_name A full path to the base model name.
#' For example, if the model with the \code{hyperparameter} block is named
#' \code{/home/user/my_model.stan}, set \code{model_name} should be equal
#' to \code{"/home/user/my_model"}
#' @param stan_data The Stan data list (e.g., for the \code{data} argument of
#' \code{stan::sampling})
#' @return A list of Stan models and parameters names which can be passed to
#' the \code{stan_sensitivity_list} argument of
#' functions like code{GetStanSensitivityFromModelFit}.
#' @export
GetStanSensitivityModel <- function(model_name, stan_data) {
    # Use a "model" with no model block to get a valid parameter list with
    # get_inits.  This way there is no worry about invalid intializations.
    model_sens_params <-
        stan(paste(model_name, "_sensitivity_parameters.stan", sep=""),
                   data=stan_data, algorithm="Fixed_param",
                   iter=1, chains=1)
    model_params <-
        stan(paste(model_name, "_generated.stan", sep=""),
                   data=stan_data, algorithm="Fixed_param",
                   iter=1, chains=1)

    sens_par_list <- get_inits(model_sens_params)[[1]]
    model_par_list <- get_inits(model_params)[[1]]

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
        model_params@.MISC$stan_fit_instance$unconstrained_param_names(
            FALSE, FALSE)
    sens_param_names <-
        model_sens_fit@.MISC$stan_fit_instance$unconstrained_param_names(
            FALSE, FALSE)

    return(list(model_sens_fit=model_sens_fit,
                model_params=model_params,
                model_sens_params=model_sens_params,
                param_names=param_names,
                sens_param_names=sens_param_names,
                sens_par_list=sens_par_list))
}


# Evaluate at draws using the hyperparameters in sens_par_list.
EvaluateAtDraws <- function(
    sampling_result, stan_sensitivity_list, sens_par_list,
    compute_grads=FALSE) {

  num_warmup_samples <- sampling_result@sim$warmup
  num_samples <- sampling_result@sim$iter - num_warmup_samples
  model_sens_fit <- stan_sensitivity_list$model_sens_fit
  sens_param_names <- stan_sensitivity_list$sens_param_names

  # Get the model gradients with respect to the hyperparameters (and parameters).
  lp_vec <- rep(NA, num_samples)
  if (compute_grads) {
    grad_mat <- matrix(NA, length(sens_param_names), num_samples)
  } else {
    grad_mat <- matrix()
  }

  cat("Evaluating sensitivity model at the MCMC draws.\n")
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


#' Process the results of Stan samples and GetStanSensitivityModel into a
#' sensitivity matrix.
#'
#' @param sampling_result The output of \code{stan::sampling}
#' @param stan_sensitivity_list The output of \code{GetStanSensitivityModel}
#' @return A list of matrices.  The elements of the list are
#' \itemize{
#'     \item{sens_mat: }{The local sensitivity of each posterior parameter to each hyperparameter.}
#'     \item{sens_mat_normalized: }{The same quantities as \code{sens_mat}, but normalized by the posterior standard deviation.}
#'     \item{grad_mat: }{The gradients of the log posterior evaluatated at the draws.}
#'     \item{lp_vec: }{The log probability of the model at each draw.}
#'     \item{draws_mat: }{The parameter draws in the same order as that of \code{grad_mat}}.
#' }
#' The result can be used directly, or passed to \code{GetTidyResult}.
#' @export
GetStanSensitivityFromModelFit <- function(
    sampling_result, stan_sensitivity_list) {

    draws_mat <- extract(sampling_result, permute=FALSE)[,1,]

    # Get the model gradients with respect to the hyperparameters (and parameters).
    model_at_draws <- EvaluateAtDraws(
            sampling_result, stan_sensitivity_list,
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
                grad_mat=grad_mat, lp_vec=model_at_draws$lp_vec,
                draws_mat=draws_mat))
}


#' Use importance sampling to approximate means at different hyperparameter
#' values.  Experimental use only.
#' @export
GetImportanceSamplingFromModelFit <- function(
    sampling_result, draws_mat, stan_sensitivity_list,
    imp_sens_par_list, lp_vec=NULL) {

    if (is.null(lp_vec)) {
        lp_vec <- EvaluateAtDraws(
            sampling_result, stan_sensitivity_list,
            stan_sensitivity_list$sens_par_list, compute_grads=FALSE)$lp_vec
    }
    # Get the model gradients with respect to the hyperparameters (and parameters).
    imp_lp_vec <- EvaluateAtDraws(
            sampling_result, stan_sensitivity_list,
            imp_sens_par_list, compute_grads=FALSE)$lp_vec

    imp_weights <- exp(imp_lp_vec - lp_vec)
    imp_weights <- imp_weights / sum(imp_weights)

    eff_num_imp_samples <- 1 / sum(imp_weights ^ 2)

    imp_means <- colSums(imp_weights * draws_mat)

    return(list(imp_weights=imp_weights, imp_lp_vec=imp_lp_vec, lp_vec=lp_vec,
                eff_num_imp_samples=eff_num_imp_samples,
                imp_means=imp_means))
}
