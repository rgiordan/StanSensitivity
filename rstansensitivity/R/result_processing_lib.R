library(ggplot2)
library(dplyr)
library(reshape2)


GetSensitivityFromGrads <- function(grad_mat, draws_mat) {
  # This should in fact match cov() but without having to transpose,
  # which gives a speedup.

  # Should match:
  # sens_mat <- cov(t(grad_mat), draws_mat)
  if (ncol(grad_mat) != nrow(draws_mat)) {
      stop(paste0(
          "The dimensions of `grad_mat` must be ",
          "# hyper parameters x # MCMC draws, and the dimensions ",
          "of draws_mat must be # draws x # parameters."))
  }

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
