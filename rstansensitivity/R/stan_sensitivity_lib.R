# Tools for getting `grad_mat` (the matrix of partial derivatives
# of the log posterior with respect to the hyperparameters) out of Stan.

library(rstan)


# Just a more readable shortcut for the Stan attribute.
GetParamNames <- function(model_fit) {
    model_fit@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)
}

#' Get the filename of the stan model to be used for sampling.
#'
#' @param model_name A full path to the base model name.
#' For example, if the model with the \code{hyperparameter} block is named
#' \code{/home/user/my_model.stan}, set \code{model_name} should be equal
#' to \code{"/home/user/my_model"}
#'
#' @return The full path of the generated file to be used for sampling.
#' @export
#' @example inst/examples/GenerateSensitivityFromModel.R
GetSamplingModelFilename <- function(model_name) {
    return(paste(model_name, "_generated.stan", sep=""))
}

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
#' @example inst/examples/GenerateSensitivityFromModel.R
GenerateSensitivityFromModel <- function(
        base_model_name,
        python_script=system.file("generate_models.py",
                                  package="rstansensitivity")) {

    model_suffix <-
        substr(base_model_name,
               nchar(base_model_name) - 4,
               nchar(base_model_name))
    stopifnot(model_suffix == ".stan")
    system(paste("python ", python_script,
                 " --base_model=", base_model_name, sep=""))
    model_name <- sub("\\.stan$", "", base_model_name)
    return(model_name)
}


#' Get a legal version of the sensitivity parameters in a list that can
#' be passed to a sensitivity model fit object.
#'
#' @param model_sens_params A stanfit object with the sensitivity model
#' parameters.  In order to guarantee legal initial values, it may be
#' preferable to use an empty model block.
#' @param model_params A stanfit object for the original model.
#' @param stan_data The stan data list for the original model with
#' hyperparameters specified for reading in the data block.  Each hyperparameter
#' in the \code{model_sens_params} stanfit model must be specified in
#' \code{stan_data}.
#' @return A list of valid model parameters than can be passed to the
#' sensivity model, in which the sampled parameters are taken from
#' \code{model_params} and the hyperparameters are taken from \code{stan_data}.
##' @export
SetSensitivityParameterList <- function(
    model_sens_params, model_params, stan_data) {

    sens_par_list <- get_inits(model_sens_params)[[1]]
    model_par_list <- get_inits(model_params)[[1]]

    for (par in names(model_par_list)) {
      sens_par_list[[par]] <- model_par_list[[par]]
    }
    for (par in setdiff(names(sens_par_list), names(model_par_list))) {
        if (!(par %in% names(stan_data))) {
            stop(sprintf("Hyperparameter %s not found in the stan_data.", par))
        }
        sens_par_list[[par]] <- stan_data[[par]]
    }
    return(sens_par_list)
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

    # This gets a model fit object for the sampling model.
    model_params <-
        stan(GetSamplingModelFilename(model_name),
             data=stan_data, algorithm="Fixed_param",
             iter=1, chains=1)

    sens_par_list <- SetSensitivityParameterList(
        model_sens_params, model_params, stan_data)

    # This is the stan fit object that we will use to evaluate the gradient
    # of the log probability at the MCMC draws.
    model_sens_fit <- rstan::stan(
        paste(model_name, "_sensitivity.stan", sep=""),
        data=stan_data, algorithm="Fixed_param",
        iter=1, chains=1, init=list(sens_par_list))

    # These names help sort through the vectors of sensitivity.
    param_names <- GetParamNames(model_params)
    sens_param_names <- GetParamNames(model_sens_fit)

    return(list(model_sens_fit=model_sens_fit,
                model_params=model_params,
                model_sens_params=model_sens_params,
                param_names=param_names,
                sens_param_names=sens_param_names,
                sens_par_list=sens_par_list))
}


#' Evaluate a model (represented as a stanfit object) at MCMC draws, possibly
#' from a different model.  The model \code{model_stanfit} is evaluated
#' at the draws found in \code{samples_stanfit}.  If there are parameters in
#' \code{model_stanfit} that are not contained in \code{samples_stanfit},
#' values from \code{model_par_list} are used.
#'
#' This was designed to evaluate gradients with respect to hyperparameters
#' which are Stan parameters in \code{model_stanfit} but not in
#' \code{samples_stanfit}.  The values of those hyperparameters would be
#' specified in \code{model_par_list} and reused for each evaluation of
#' \code{model_stanfit}.  Since this is the intended use case, the function
#' requires that \code{model_stanfit} contains a (non-strict) superset of the
#' parameters in \code{samples_stanfit}.
#'
#' @param samples_stanfit (sampling_result) The output of Stan sampling
#' containing the draws at which you want to evaluate the model.
#' (e.g. \code{rstan::sampling()})
#' @param model_stanfit (model_fit) A stanfit object from a model, possibly containing a
#' superset of the parameters in samples_stanfit.
#' @param model_par_list A list of parameters for model_stanfit which can be
#' passed to
#' \code{unconstrain_pars} for the \code{model_stanfit} model.  Every parameter
#' in \code{get_inits} applied to \code{samples_stanfit} must be also be found
#' in \code{model_par_list}.
#' @param compute_grads If FALSE, only return the log probability.
#' @param max_chains How many chains to evaluate.  The default is to
#' evaluate the number of chains in \code{samples_stanfit}.
#' @param max_num_samples How many samples to evaluate.  The default is to
#' evaluate the number of non-warmup samples in \code{samples_stanfit}.
#' @param verbose If true, display progress messages.
#' @return
#' A list with \code{lp_vec} containing a vector of log probabilities and,
#' if \code{compute_grads} is \code{TRUE}, gradients of the log probability,
#' all with resepct to the parameters in \code{model_stanfit} at the draws in
#' \code{samples_stanfit}.  If the sampler contains multiple chains the
#' chains are concatenated in order.
##' @export
EvaluateModelAtDraws <- function(
    samples_stanfit, model_stanfit, model_par_list,
    compute_grads=FALSE, max_chains=Inf, max_num_samples=Inf,
    verbose=TRUE) {

  num_warmup_samples <- samples_stanfit@sim$warmup
  num_samples <- min(samples_stanfit@sim$iter - num_warmup_samples,
                     max_num_samples)
  num_chains <- min(samples_stanfit@sim$chains, max_chains)

  # Check that every parameter in the samples_stanfit is also in model_stanfit.
  samples_stanfit_pars <- rstan::get_inits(samples_stanfit, iter=1)[[1]]
  missing_pars <- setdiff(names(samples_stanfit_pars), names(model_par_list))
  if (length(missing_pars) > 0) {
      err_msg <- paste0(
          "Every parameter in `samples_stanfit` must also be found in ",
          "`model_stanfit`.  ",
          sprintf(
              "Parameters from sampling result not in new model: %s",
              paste(missing_pars, collapse=", ")))
      stop(err_msg)
  }

  # Get the model gradients with respect to the hyperparameters (and parameters).
  lp_vec <- rep(NA, num_samples)
  if (compute_grads) {
    param_names <- GetParamNames(model_stanfit)
    grad_mat <- matrix(
        NA, nrow=length(param_names), ncol=num_samples * num_chains)
    rownames(grad_mat) <- param_names
  } else {
    grad_mat <- matrix()
  }

  for (chain in 1:num_chains) {
      if (verbose) {
          cat("Evaluating model at the MCMC draws for chain ", chain, ".\n")
          prog_bar <- txtProgressBar(min=1, max=num_samples, style=3)
      }
      for (n in 1:num_samples) {
         if (verbose) {
             setTxtProgressBar(prog_bar, value=n)
         }

        # We rely on get_inits to return the draws at iteration n in a form
        # that is easy to parse.
        par_list <- get_inits(
            samples_stanfit, iter=n + num_warmup_samples)[[chain]]
        for (par in ls(par_list)) {
          model_par_list[[par]] <- par_list[[par]]
        }
        pars_free <- rstan::unconstrain_pars(model_stanfit, model_par_list)

        # The index in the matrix stacked by chain.
        # This needs to match the stacking done to the Stan samples in
        # StackChainArray().  Perhaps it would be wiser to store the gradients
        # in the same order as rstan::extract() and use the same function to
        # put them into compatible shapes...
        ix <- (chain - 1) * num_samples + n
        if (compute_grads) {
          glp <- rstan::grad_log_prob(model_stanfit, pars_free)
          grad_mat[, ix] <- glp
          lp_vec[ix] <- attr(glp, "log_prob")
        } else {
          lp_vec[ix] <- rstan::log_prob(model_stanfit, pars_free)
        }
      }
      if (verbose) {
          close(prog_bar)
      }
  }
  return(list(lp_vec=lp_vec, grad_mat=grad_mat))
}


#' Check whether two stan models objects agree on draws in
#' \code{sampling_result}.
#'
#' @param sampling_result Samples drawn from a Stan model, possibly distinct
#' from \code{model_1} and \code{model_2}.
#' @param model_1 The path to the first stan model.
#' @param model_2 The path to the second stan model.
#' @param stan_data_1 A stan data file for the first model.
#' @param stan_data_2 A stan data file for the second model.
#' @param check_grads Whether to check the gradients of the log probabilities.
#' @param tol The tolerance for the difference in grads and log probability.
#' @return A boolean indicating whether model_1 and model_2 agree with each
#' other on the samples specified in sampling_result.  Log probability is
#' is compared up to a constant, and the gradients are compared on the
#' parameters in common according to the unconstrained parameter names.
#' @export
CheckModelEquivalence <- function(
    sampling_result, model_1, model_2, stan_data_1, stan_data_2,
    check_grads=TRUE,
    tol=1e-8, max_chains=1, max_num_samples=100) {

    model_fit_1 <- rstan::sampling(
        object=model_1, data=stan_data_1,
        algorithm="Fixed_param", iter=1, chains=1)
    model_par_list_1 <- get_inits(model_fit_1, 1)[[1]]

    model_fit_2 <- rstan::sampling(
        object=model_2, data=stan_data_2,
        algorithm="Fixed_param", iter=1, chains=1)
    model_par_list_2 <- get_inits(model_fit_2, 1)[[1]]

    return(CheckModelFitEquivalence(
                sampling_result, model_fit_1, model_fit_2,
                model_par_list_1, model_par_list_2,
                check_grads=check_grads, tol=tol,
                max_chains=max_chains, max_num_samples=max_num_samples))
}


CheckModelFitEquivalence <- function(
    sampling_result,
    model_fit_1, model_fit_2,
    model_par_list_1, model_par_list_2,
    check_grads=TRUE,
    tol=1e-8, max_chains=1, max_num_samples=100) {

    model_draws_1 <- EvaluateModelAtDraws(
        sampling_result, model_fit_1, model_par_list_1,
        compute_grads=check_grads,
        max_chains=max_chains, max_num_samples=max_num_samples)
    model_draws_2 <- EvaluateModelAtDraws(
        sampling_result, model_fit_2, model_par_list_2,
        compute_grads=check_grads,
        max_chains=max_chains, max_num_samples=max_num_samples)

    # The log probability can differ up to a constant.
    log_prob_ok <-
        sqrt(var(model_draws_1$lp_vec - model_draws_2$lp_vec)) < tol

    if (check_grads) {
        # Ideally we would compare on the unconstrained parameter names used
        # to generate the sampling result.  Unfortunately, it appears that if
        # a sampling result has been saved to an Rdata file,
        # the information necessary to run GetParamNames is not saved.
        common_par_names <-
            intersect(rownames(model_draws_1$grad_mat),
                      rownames(model_draws_2$grad_mat))
        stopifnot(length(common_par_names) > 0)
        grad_mat_1 <- model_draws_1$grad_mat[common_par_names, ]
        grad_mat_2 <- model_draws_2$grad_mat[common_par_names, ]
        grad_ok <- max(abs(grad_mat_1 - grad_mat_2)) < tol
    } else {
        grad_ok <- TRUE
    }

    return(log_prob_ok & grad_ok)
}


# Evaluate at draws using the hyperparameters in sens_par_list.
EvaluateAtDraws <- function(
    sampling_result, stan_sensitivity_list, sens_par_list,
    compute_grads=FALSE) {

    return(EvaluateModelAtDraws(
        sampling_result,
        stan_sensitivity_list$model_sens_fit,
        stan_sensitivity_list$sens_par_list,
        compute_grads=compute_grads))
}


# extract() returns an array of samples x chain x parameter.  Stack them
# by chain into a single matrix of (samples * chains) x parameter.
StackChainArray <- function(draws_array) {
    num_chains <- dim(draws_array)[2]
    draws_mat <- do.call(
        rbind, lapply(1:num_chains, function(chain) { draws_array[, chain, ] }))
    return(draws_mat)
}


#' Process the results of Stan samples and GetStanSensitivityModel into a
#' sensitivity matrix.  Note that currently, only the first chain is supported.
#'
#' @param sampling_result The output of \code{stan::sampling}
#' @param stan_sensitivity_list The output of \code{GetStanSensitivityModel}
#' @return A list of matrices.  The elements of the list are
#' \itemize{
#'     \item{grad_mat: }{The gradients of the log posterior evaluatated at the draws.}
#'     \item{lp_vec: }{The log probability of the model at each draw.}
#'     \item{draws_mat: }{The parameter draws in the same order as that of \code{grad_mat}}.
#' }
#' @export
GetStanSensitivityMatricesFromModelFit <- function(
    sampling_result, stan_sensitivity_list) {

    draws_mat <- StackChainArray(extract(sampling_result, permute=FALSE))

    # Get the model gradients with respect to the hyperparameters (and parameters).
    model_at_draws <- EvaluateAtDraws(
            sampling_result, stan_sensitivity_list,
            stan_sensitivity_list$sens_par_list, compute_grads=TRUE)

    # Stan takes gradients with respect to everything in the parameters block,
    # not just the hyperparameters.  Remove the rows not corresponding to
    # hyperparameters.
    param_names <- stan_sensitivity_list$param_names
    sens_param_names <- stan_sensitivity_list$sens_param_names
    grad_mat <- model_at_draws$grad_mat[
        setdiff(sens_param_names, param_names),, drop=FALSE]

    return(list(grad_mat=grad_mat,
                lp_vec=model_at_draws$lp_vec,
                draws_mat=draws_mat))
}


#' Process the results of Stan samples and GetStanSensitivityModel into a
#' sensitivity matrix.  Note that currently, only the first chain is supported.
#'
#' @param sampling_result The output of \code{stan::sampling}
#' @param stan_sensitivity_list The output of \code{GetStanSensitivityModel}
#' @return A list of matrices.  The elements of the list are
#' \itemize{
#'     \item{sens_mat: }{The local sensitivity of each posterior parameter to each hyperparameter.}
#'     \item{sens_mat_normalized: }{\code{sens_mat}, but normalized by the posterior standard deviation.}
#'     \item{grad_mat: }{The gradients of the log posterior evaluatated at the draws.}
#'     \item{lp_vec: }{The log probability of the model at each draw.}
#'     \item{draws_mat: }{The parameter draws in the same order as that of \code{grad_mat}}.
#' }
#' The result can be used directly, or passed to \code{GetTidyResult}.
#' @export
GetStanSensitivityFromModelFit <- function(
    sampling_result, stan_sensitivity_list) {

    # TODO: deprecate this for something more general.
    sens_mats <- GetStanSensitivityMatricesFromModelFit(
        sampling_result, stan_sensitivity_list)

    # Calculate the sensitivity.
    sens_mats$sens_mat <- GetSensitivityFromGrads(
        sens_mats$grad_mat, sens_mats$draws_mat)

    # Normalize by the marginal standard deviation.
    sens_mats$sens_mat_normalized <- NormalizeSensitivityMatrix(
        sens_mats$sens_mat, sens_mats$draws_mat)

    return(sens_mats)
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

    return(list(imp_weights=imp_weights,
                imp_lp_vec=imp_lp_vec,
                lp_vec=lp_vec,
                eff_num_imp_samples=eff_num_imp_samples,
                imp_means=imp_means))
}
