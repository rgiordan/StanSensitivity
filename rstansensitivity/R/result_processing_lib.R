library(ggplot2)
library(dplyr)
library(reshape2)


#' Convert a covariance matrix to a dataframe.
#'
#' @param cov_mat The covariance matrix.
#' @param remove_repeats Optional.  If TRUE, and if cov_mat
#' is symmetric, include only the upper triangular part of cov_mat.
#' Default value = FALSE.
#' @return A data frame with columns row_variable, column_variable, and
#' value, where value is the specificed covariance.
#' @export
CovarianceMatrixToDataframe <- function(cov_mat, remove_repeats=FALSE) {
  if (is.null(colnames(cov_mat))) {
    column_names <- paste0("col", 1:ncol(cov_mat))
    colnames(cov_mat) <- column_names
  }
  if (is.null(rownames(cov_mat))) {
    row_names <- paste0("row", 1:nrow(cov_mat))
    rownames(cov_mat) <- row_names
  }
  cov_df <- data.frame(cov_mat, stringsAsFactors=FALSE)
  names(cov_df) <- column_names
  row_df <- data.frame(row_variable=row_names, row=1:length(row_names),
                       stringsAsFactors=FALSE)
  col_df <- data.frame(column_variable=column_names, col=1:length(column_names),
                       stringsAsFactors=FALSE)
  cov_df <- cov_df %>%
           mutate(row_variable=row_names) %>%
           melt(id.var="row_variable") %>%
           rename(column_variable=variable)
  if (remove_repeats) {
    if (length(row_names) != length(column_names)) {
      stop(paste0("To use remove_repeats, the row names and column names ",
                  "must be the same length."))
    }
    if (any(row_names != column_names)) {
      stop(paste0("To use remove_repeats, the row names and column names must ",
                  "be identical"))
    }
    if (max(abs(cov_mat - t(cov_mat)), na.rm=TRUE) > 1e-8) {
      stop("To use remove_repeats, cov_mat must be symmetric.")
    }
    cov_df <-
      cov_df %>%
      inner_join(row_df, by="row_variable") %>%
      inner_join(col_df, by="column_variable") %>%
      filter(row <= col) %>%
      select(-row, -col)
  }
  return(cov_df)
}



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


#' Make a tidy dataframe out of a sensitivity matrix and its standard errors.
#' @param sens_result The output of \code{GetStanSensitivityFromModelFit}.
#' @param transform_matrix A matrix where each row is a linear combination
#' of parameters for which to calculate sensitivity.  The row names should
#' describe the linear combinations.
#' @return A new sensitivity result list for the specified linear combinations.
#' @export
TransformSensitivityResult <- function(sens_result, transform_matrix) {
    stopifnot(ncol(transform_matrix) == nrow(sens_result$sens_mat))
    if (is.null(rownames(transform_matrix))) {
        rownames(transform_matrix) <-
            paste("transform", 1:nrow(transform_matrix), sep="_")
    }
    transform_sens_result <- sens_result
    transform_sens_result$sens_mat <-
      transform_matrix %*% transform_sens_result$sens_mat
    # The columns are normalized so this preserves the normalization.
    transform_sens_result$sens_mat_normalized <-
      transform_matrix %*% transform_sens_result$sens_mat_normalized
    transform_sens_result$grad_mat <-
      transform_matrix %*% transform_sens_result$grad_mat
    return(transform_sens_result)
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
#' @param sens_mat A matrix of sensitivities.
#' @param sens_se A matrix of standard errors of \code{sens_mat}.
#' @param measure What to call these sensitivites.
#' @param num_se The number of standard errors for the upper and lower bounds.
#' @return A tidy dataframe with columns for the parameters, hyperparameters,
#' sensitivities, and their standard errors.
#' @export
SummarizeSensitivityMatrices <- function(
    sens_mat, sens_se, measure, num_se=2) {

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
        sens_df,sens_norm_df, by=c("hyperparameter", "parameter"))

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
