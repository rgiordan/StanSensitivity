draws1_matdraws1_matdraws1_matlibrary(mcmcse)


mcse.multi_safe <- function(arg_draws) {
    # Set a default to return if the method fails.
    output_list <- list(cov=matrix(NA, ncol(arg_draws), ncol(arg_draws)))
    tryCatch(output_list <- mcse.multi_safe(arg_draws),
             error=function(e) { print(e) },
             warning=function(w) { print(w) }) # Should we suppress the warning?
    return(output_list$cov)
}


GetCovarianceSE <- function(x_draws, y_draws, fix_mean=FALSE) {
  x_mean <- mean(x_draws)
  y_mean <- mean(y_draws)
  if (fix_mean) {
    g_draws <- x_draws * y_draws - x_mean * y_mean
    g_se <- mcmcse::mcse(g_draws)$se
  } else {
    arg_draws <- cbind(x_draws * y_draws, x_draws, y_draws)
    arg_cov_mat <- mcse.multi_safe(arg_draws)$cov
    grad_g <- c(1, -1 * y_mean, -1 * x_mean)
    g_se <- as.numeric(sqrt(t(grad_g) %*% arg_cov_mat %*% grad_g / nrow(arg_draws)))
  }
  return(g_se)
}


PackNormalizedCovariancePar <- function(x_draws, y_draws) {
  return(cbind(
    x_draws * y_draws,
    x_draws ^ 2,
    x_draws,
    y_draws
  ))
}


GetNormalizedCovarianceGradient <- function(par) {
  mean_xy <- par[1]
  mean_x2 <- par[2]
  mean_x <- par[3]
  mean_y <- par[4]
  sigma_xy <- mean_xy - mean_x * mean_y
  sigma_x <- sqrt(mean_x2 - mean_x^2)
  return(c(
    1 / sigma_x,
    -0.5 * sigma_xy / (sigma_x^3),
    (sigma_xy * mean_x / (sigma_x^2) - mean_y) / sigma_x,
    -mean_x / sigma_x
  ))
}


# For testing.
GetNormalizedCovariance <- function(par) {
  mean_xy <- par[1]
  mean_x2 <- par[2]
  mean_x <- par[3]
  mean_y <- par[4]
  sigma_xy <- mean_xy - mean_x * mean_y
  sigma_x <- sqrt(mean_x2 - mean_x^2)
  return(sigma_xy / sigma_x)
}


GetNormalizedCovarianceSE <- function(x_draws, y_draws) {
  par_draws <- PackNormalizedCovariancePar(x_draws, y_draws)
  par_cov_mat <- mcse.multi_safe(par_draws)$cov
  grad_g <- GetNormalizedCovarianceGradient(colMeans(par_draws))
  return(
    as.numeric(sqrt(t(grad_g) %*% par_cov_mat %*% grad_g / nrow(par_draws))))
}


#' Estimate Monte Carlo standard errors for the sensitivity matrices using
#' a normal approximation and the delta method.  These standard errors are
#' based on what may be unrealistic assumptions, and should be interpreted with
#' caution.
#'
#' @param draws_mat Parameter draws from \code{GetStanSensitivityFromModelFit}.
#' @param grad_mat Gradients from \code{GetStanSensitivityFromModelFit}.
#' @param fix_mean Experimental: estimate standard errors as if the means were
#' known exactly.  Not recommended.
#' @param normalized Compute standard errors for the normalized sensitivities.
#' @return A matrix of the same dimensions as \code{sens_mat} containing
#' estimates of the Monte Carlo standard errors of the sensitivities.
#' @export
GetSensitivityStandardErrors <- function(
    draws_mat, grad_mat, fix_mean=FALSE, normalized=FALSE) {

  if (normalized && fix_mean) {
    stop("You cannot specify both fix_mean and normalized.")
  }
  sens_param_names <- rownames(grad_mat)
  param_names <- colnames(draws_mat)

  cov_se_mat <- matrix(NA, nrow(grad_mat), ncol(draws_mat))
  rownames(cov_se_mat) <- sens_param_names
  colnames(cov_se_mat) <- param_names
  for (par_ind in 1:length(param_names)) {
    for (grad_ind in 1:length(sens_param_names)) {
      par_draws <- draws_mat[, param_names[par_ind]]
      grad_draws <- grad_mat[sens_param_names[grad_ind], ]
      if (normalized) {
        cov_se_mat[grad_ind, par_ind] <-
          GetNormalizedCovarianceSE(par_draws, grad_draws)
      } else {
        cov_se_mat[grad_ind, par_ind] <-
          GetCovarianceSE(par_draws, grad_draws, fix_mean=fix_mean)
      }
    }
  }
  return(cov_se_mat)
}



# Bootstrap the covariance matrices using a block bootstrap on the
# MCMC chain.
GetBlockBootstrapStandardErrors <- function(draws1_mat, draws2_mat,
                                            num_blocks, num_draws) {

    if (nrow(draws1_mat) != nrow(draws2_mat)) {
        stop("draws1_mat and draws2_mat must have the same number of rows.")
    }

    num_samples <- nrow(draws1_mat)

    block_size <- floor(num_samples / num_blocks)

    # Correction factor if the number of blocked observations is not the same
    # as the original.
    n_factor <- (block_size * num_blocks) / num_samples

    # The indices of each block into the MCMC samples.
    block_inds <- lapply(
      1:num_blocks,
      function(ind) { (ind - 1) * block_size + 1:block_size })

    # Lists of the draws of each block.
    draws1_blocks <- lapply(
      1:num_blocks, function(ind) { draws1_mat[block_inds[[ind]], ]} )
    draws2_blocks <- lapply(
      1:num_blocks, function(ind) { draws2_mat[block_inds[[ind]], ]} )

    # Bind blocks together according to random block indices.
    DrawBlocks <- function(draw_blocks, block_ind_draws) {
        return(do.call(rbind,
                       lapply(block_ind_draws,
                              function(ind) { draw_blocks[[ind]] })))
    }

    cov_samples <- array(NA, c(num_draws, ncol(draws1_mat), ncol(draws2_mat)))
    for (draw in 1:num_draws) {
        block_ind_draws <- sample(1:num_blocks, num_blocks, replace=TRUE)

        draws1_mat_sampled <- DrawBlocks(draws1_blocks, block_ind_draws)
        draws2_mat_sampled <- DrawBlocks(draws2_blocks, block_ind_draws)

        cov_samples[draw, , ] <-
          ComputeIJCovariance(draws1_mat_sampled, draws2_mat_sampled)
    }

    cov_se <- sqrt(n_factor) * apply(cov_samples, MARGIN=c(2, 3), sd)

    return(list(cov_samples=cov_samples, cov_se=cov_se))
}
