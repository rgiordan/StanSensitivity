# Functions to process stan and sensitivity results into tidy dataframes.
#
# TODO: move most of this over to tidybayes

library(ggplot2)
library(dplyr)
library(reshape2)


#' Get a dataframe with the hyperparameter values contained in a stan data list.
#'
#' @param stan_sensitivity_list The output of \code{GetStanSensitivityModel}
#' @param stan_data The Stan data for the original model.
#' @return A tidy dataframe with the hyperparameter names and values.
#' @export
GetHyperparameterDataFrame <- function(stan_sensitivity_list, stan_data) {
  sens_par_list <- SetSensitivityParameterList(
    stan_sensitivity_list$model_sens_params,
    stan_sensitivity_list$model_params,
    stan_data)
  pars_free <- rstan::unconstrain_pars(
      stan_sensitivity_list$model_sens_fit, sens_par_list)
  names(pars_free) <- stan_sensitivity_list$sens_param_names
  hyperparam_names <- setdiff(stan_sensitivity_list$sens_param_names,
                              stan_sensitivity_list$param_names)
  hyperparameter_df <-
    data.frame(hyperparameter=hyperparam_names,
               hyperparameter_val=pars_free[hyperparam_names])
  rownames(hyperparameter_df) <- NULL
  return(hyperparameter_df)
}


#' Execute an R file and store its environment variables in a list.
#' @param data_file The filename of a code containing R code..
#' @return A list with the variables defined when running `data_file`.
#' @export
LoadStanData <- function(data_file) {
  stan_data <- new.env()
  source(data_file, local=stan_data)
  stan_data <- as.list(stan_data)
  return(stan_data)
}

SensitivityMatrixToDataframe <- function(
    sens_mat, hyperparameter_names, parameter_names) {
  colnames(sens_mat) <- parameter_names
  return(data.frame(sens_mat) %>%
    mutate(hyperparameter=hyperparameter_names) %>%
    melt(id.var="hyperparameter") %>%
    rename(parameter=variable))
}


#' Make a tidy dataframe out of a sensitivity matrix and its standard errors.
#' @param sens_mat A matrix of sensitivities.  Hyperparameters should be in the rows and parameters in the columns.
#' @param sens_se A matrix of standard errors of \code{sens_mat}.
#' @param measure What to call these sensitivites.
#' @param num_se The number of standard errors for the upper and lower bounds.
#' @return A tidy dataframe with columns for the parameters, hyperparameters,
#' sensitivities, and their standard errors.
#' @export
SummarizeSensitivityMatrices <- function(
    sens_mat, sens_se, measure, num_se=qnorm(0.975)) {

    sens_df <- rbind(
      SensitivityMatrixToDataframe(
        sens_mat,
        hyperparameter_names=rownames(sens_mat),
        parameter_names=colnames(sens_mat)) %>%
        mutate(measure=measure),
      SensitivityMatrixToDataframe(
        sens_se,
        hyperparameter_names=rownames(sens_se),
        parameter_names=colnames(sens_se)) %>%
        mutate(measure=paste(measure, "se", sep="_")),
      SensitivityMatrixToDataframe(
        sens_mat + num_se * sens_se,
        hyperparameter_names=rownames(sens_mat),
        parameter_names=colnames(sens_mat)) %>%
        mutate(measure=paste(measure, "upper", sep="_")),
      SensitivityMatrixToDataframe(
        sens_mat - num_se * sens_se,
        hyperparameter_names=rownames(sens_mat),
        parameter_names=colnames(sens_mat)) %>%
        mutate(measure=paste(measure, "lower", sep="_"))
      ) %>%
      dcast(hyperparameter + parameter ~ measure) %>%
      filter(parameter != "lp__")
    return(sens_df)
}

#' Process the results of GetStanSensitivityFromModelFit into a
#' tidy dataframe with standard errors.
#'
#' @param sens_result The output of \code{GetStanSensitivityFromModelFit}.
#' @param num_se The number of standard errors for the upper and lower bounds.
#' @return A dataframe summarizing the sensitivity of the model posterior means
#' to the hyperparameters.  The reported standard errors are based on a
#' multivariate normal and delta method approximation which may not be
#' accurate for highly variable or non-normal draws. The standard errors should
#' be interpreted with caution.
#' \itemize{
#'     \item{hyperparameter: }{The name of the hyperparameter.}
#'     \item{parameter: }{The name of the parameter.}
#'     \item{sensitivity:}
#'     {The MCMC estimate of \code{dE[parameter | X, hyperparameter] /  d hyperparameter}.}
#'     \item{sensitivity_se: }
#'     {The estimated Monte Carlo error of \code{sensitivity}.}
#'     \item{normalized_sensitivity: }
#'     {The MCMC estimate of \code{sensitivity / sd(parameter}).}
#'     \item{normalized_sensitivity_se: }
#'     {The estimated Monte Carlo error of \code{normalized_sensitivity}.}
#' }
#' @export
GetTidyResult <- function(sens_result, num_se=2) {
    sens_se <- GetSensitivityStandardErrors(
        sens_result$draws_mat, sens_result$grad_mat,
        fix_mean=FALSE, normalized=FALSE)
    norm_sens_se <- GetSensitivityStandardErrors(
        sens_result$draws_mat, sens_result$grad_mat,
        fix_mean=FALSE, normalized=TRUE)
    sens_df <- SummarizeSensitivityMatrices(
        sens_result$sens_mat, sens_se,
        measure="sensitivity", num_se=num_se)
    sens_norm_df <- SummarizeSensitivityMatrices(
        sens_result$sens_mat_normalized, norm_sens_se,
        measure="normalized_sensitivity", num_se=num_se)

    result <- inner_join(
        sens_df, sens_norm_df, by=c("hyperparameter", "parameter"))

  return(result)
}


#' Make a barchart from the output of \code{GetTidyResult}.
#'
#' @param sens_df The output of \code{GetTidyResult}.
#' @param normalized Whether to display the normalized sensitivities.
#' @param se_num The number of standard errors to plot with the errorbars.
#' @return A \code{ggplot} object showing the sensitivites of each parameter
#' to each hyperparameter.
#' @export
PlotSensitivities <- function(sens_df, normalized=TRUE, se_num=2) {
    if (normalized) {
        y_axis_label <- "Normalized sensitivity"
        sens_df_plot <-
            mutate(sens_df, s=normalized_sensitivity,
                   s_se=normalized_sensitivity_se)
    } else {
        y_axis_label <- "Unnormalized sensitivity"
        sens_df_plot <-
            mutate(sens_df, s=sensitivity, s_se=sensitivity_se)
    }
    return(
    ggplot(sens_df_plot) +
      geom_bar(aes(x=parameter, y=s, fill=hyperparameter),
                   stat="identity", position="dodge") +
      geom_errorbar(aes(x=parameter,
                        ymin=s - se_num  * s_se,
                        ymax=s + se_num  * s_se,
                        group=hyperparameter),
                    position=position_dodge(0.9), width=0.2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_discrete(name="Hyperparameter") +
      ylab(y_axis_label) + xlab("Parameter")
    )
}


#' Convert a Stan sampling result to a form that can be joined with tidy
#' sensitivity results.
#'
#' @param sampling_result The output of \code{stan::sampling}
#' @param cols The columns of the summary to keep.
#' @return A data frame with the MCMC reuslts that can be joined with tidy
#' sensitivity results.
#' @export
GetMCMCDataFrame <- function(
    sampling_result, cols=c("mean", "se_mean", "sd", "n_eff", "Rhat")) {
  mcmc_result <- as.data.frame(rstan::summary(sampling_result)$summary)
  mcmc_result$parameter <- make.names(rownames(mcmc_result))
  rownames(mcmc_result) <- NULL
  mcmc_result <-
    select(mcmc_result, "parameter", cols) %>%
    filter(parameter != "lp__")
  return(mcmc_result)
}


#' Use the linear approximation to predict the sensitivity to a new
#' stan data list.
#'
#' @param stan_sensitivity_list The output of \code{GetStanSensitivityModel}
#' @param stan_result The output of \code{GetStanSensitivityFromModelFit}
#' @param stan_data The original stan data at which
#' \code{stan_sensitivity_list} was calculated.
#' @param stan_data_perturb A new stan data file with different hyperparameters.
#' @param description A hyperparameter name to describe this perturbation.
#' @return A tidy sensitivity dataframe where the sensitivity is in the
#' direction of the difference between the hyperparameters in the two stan
#' data lists.
#' @export
PredictSensitivityFromStanData <- function(
    stan_sensitivity_list, sens_result, stan_data, stan_data_perturb,
    description="perturbation") {

    hyperparameter_df <-
      inner_join(
          GetHyperparameterDataFrame(stan_sensitivity_list, stan_data) %>%
               rename(hyperparameter_val_orig=hyperparameter_val),
          GetHyperparameterDataFrame(stan_sensitivity_list, stan_data_perturb),
               by="hyperparameter") %>%
      mutate(hyperparameter_diff=hyperparameter_val - hyperparameter_val_orig)
    linear_comb <- matrix(hyperparameter_df$hyperparameter_diff, nrow=1)
    colnames(linear_comb) <- hyperparameter_df$hyperparameter
    linear_comb <- linear_comb[, rownames(sens_result$sens_mat), drop=FALSE]
    rownames(linear_comb) <- description
    sens_result_pert <- TransformSensitivityResult(sens_result, linear_comb)
    return(GetTidyResult(sens_result_pert))
}


# Get a tidy version of a Stan summary with tidybayes.
#
#' @param stanfit A \code{stanfit} object (e.g., from \code{sampling}).
#' @param ... Parameter names in the style of \code{tidybayes}.
#' @param spread Optional.  If \code{TRUE}, return a wide dataframe in which
#'        each summary metric is a column.  Otherwise, return a tall dataframe
#'        in which the summary metric is in a column called \code{metric}.
#' @return A tidy dataframe containing the \code{stan} summary object.
#' @export
GetTidyStanSummary <- function(stanfit, ..., spread=FALSE) {
  pars <- enquos(...)
  summary_mat <- t(rstan::summary(stanfit)$summary)
  tidy_summary <-
    gather_draws(summary_mat, !!!pars) %>%
    inner_join(data.frame(.draw=1:nrow(summary_mat),
                          metric=rownames(summary_mat)),
               by=".draw") %>%
    select(-.chain, -.draw, -.iteration)

  if (spread) {
    key_string <- paste(setdiff(names(tidy_summary),
                        c("metric", ".value")), collapse=" + ")
    tidy_summary <- dcast(
        tidy_summary,
        formula(sprintf("%s ~ metric", key_string)),
        value.var=".value")
  }
  return(tidy_summary)
}


RemoveExtraTidyColumns <- function(df) {
  select(df, -.chain, -.iteration, -.draw)
}


#' Convert a column of variable names to a tidy format.
#' @export
TidyColumn <- function(df, col, ..., variable_name=".variable") {
  pars <- enquos(...)
  par_names <- unique(df[[col]])
  tidy_pars_df <-
    matrix(par_names, nrow=1, dimnames=list("", par_names)) %>%
    gather_draws(!!!pars) %>%
    rename(!!sym(col):=.value, !!sym(variable_name):=.variable) %>%
    select(-.chain, -.iteration, -.draw)
  return(inner_join(df, tidy_pars_df, by=col))
}


GatherRowNamedMatrix <- function(mat, ...) {
    pars <- enquos(...)
    df <-
        tidybayes::gather_draws(mat, !!!pars) %>%
        inner_join(data.frame(.draw=1:nrow(mat), rowname=rownames(mat)),
            by=".draw") %>%
        RemoveExtraTidyColumns()
    return(df)
}


#' @export
GetTidySensitivityResults <- function(grad_mat,
                                      draws_mat,
                                      ...,
                                      normalize=TRUE,
                                      calculate_se=TRUE) {
  pars <- enquos(...)
  GetSensitivityDf <- function(mat, measure_name) {
    GatherRowNamedMatrix(mat, !!!pars) %>%
      mutate(measure=measure_name) %>%
      rename(hyperparameter=rowname)
  }
  sens_mat <- GetSensitivityFromGrads(grad_mat, draws_mat)
  results <- GetSensitivityDf(sens_mat, "sensitivity")
  if (calculate_se) {
    sens_mat_se <- GetSensitivityStandardErrors(
      draws_mat, grad_mat, normalized=FALSE)
    results <- bind_rows(
      results,
      GetSensitivityDf(sens_mat_se, "sensitivity_se"))
  }

  if (normalize) {
    # Normalize by the marginal standard deviation.
    sens_mat_norm <- NormalizeSensitivityMatrix(sens_mat, draws_mat)
    results <- bind_rows(
      results,
      GetSensitivityDf(sens_mat_norm, "normalized_sensitivity"))
    if (calculate_se) {
      sens_mat_norm_se <- GetSensitivityStandardErrors(
        draws_mat, grad_mat, normalized=TRUE)
      results <- bind_rows(
        results,
        GetSensitivityDf(sens_mat_norm_se, "normalized_sensitivity_se"))
    }
  }
  return(results)
}
