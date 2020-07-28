#' Aggregate the draws of the log likelihood according to a grouping variable
#' where the observations are exchangeable within a group.
#' @param loglik_draws_mat Draws of the log likelihood.
#' @param exchangeable_col Optional, a vector of indices indicating the grouping
#' of exchangeable observations.
#' @return The log likelihood draws of aggregated by the given grouping.
#' @export
GroupLogLikelihoodDraws <- function(loglik_draws_mat, exchangeable_col) {
  if (length(exchangeable_col) != ncol(loglik_draws_mat)) {
    stop(paste0("exchangeable_col must be as long as the number of ",
                "columns in loglik_draws_mat."))
  }

  # Sum the log likelihood within exchangable observations, drop the
  # grouping column, and re-cast as a matrix.
  loglik_draws_mat <- aggregate(
      t(loglik_draws_mat),
      by=list(exchangeable_col=exchangeable_col), sum)[, -1] %>%
    t() %>%
    as.matrix()
  return(loglik_draws_mat)
}



#' Return an estimate of the infinitesimal jackknife covariance estimate
#' of the frequentist variance of the posterior expectations of the parameters
#' in param_draws_mat.
#' @param loglik_draws_mat Draws of the log likelihood.
#' @param param_draws_mat Draws of the parameters of interest..
#' @return The IJ covariance matrix, which is an estimate of
#' N * Cov(E[params | x]), where N is the number of distinct exchangeable
#' observations.
#' @export
ComputeIJCovariance <- function(loglik_draws_mat, param_draws_mat,
                                exchangeable_col=NULL) {
  if (nrow(loglik_draws_mat) != nrow(param_draws_mat)) {
    stop(paste0("loglik_draws_mat and param_draws_mat must have the ",
                "same number of rows."))
  }
  num_obs <- ncol(loglik_draws_mat)
  infl_draws_mat <- num_obs * cov(loglik_draws_mat, param_draws_mat)
  colnames(infl_draws_mat) <- colnames(param_draws_mat)
  ij_cov <- cov(infl_draws_mat, infl_draws_mat)
  return(ij_cov)
}
