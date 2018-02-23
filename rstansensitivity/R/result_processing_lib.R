library(ggplot2)
library(dplyr)
library(reshape2)


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
    stopifnot(ncol(transform_matrix) == nrow(sens_mat))
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
