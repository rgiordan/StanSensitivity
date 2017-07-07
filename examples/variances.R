# Compute the variance of the covariances and their ratios.

GetSensitivityFromGrads <- function(grad_mat, draws_mat) {
  grad_means <- rowMeans(grad_mat)
  draw_means <- colMeans(draws_mat)
  
  n <- nrow(draws_mat)
  sens_mat <- (grad_mat %*% draws_mat) / (n - 1) - grad_means %*% t(draw_means) * n / (n - 1)
  return(sens_mat)
}

NormalizeSensitivityMatrix <- function(sens_mat, draws_mat) {
  draw_sds <- sqrt(apply(draws_mat, sd, MARGIN=2))
  return(sens_mat / rep(draw_sds, each=nrow(sens_mat)))
}

param_names <- stan_sensitivity_list$param_names
sens_param_names <- stan_sensitivity_list$sens_param_names
grad_mat <- sens_result$grad_mat

rownames(grad_mat) <- sens_param_names
grad_mat <- grad_mat[setdiff(sens_param_names, param_names), , drop=FALSE]

num_boot <- 200
sens_mat <- GetSensitivityFromGrads(grad_mat, draws_mat)
sens_mat_array <- array(NA, dim=c(num_boot, dim(sens_mat)))
sens_mat_normalized_array <- array(NA, dim=c(num_boot, dim(sens_mat)))

cat("Bootstrapping sensitivity matrix.")
prog_bar <- txtProgressBar(min=1, max=num_boot, style=3)
for (boot in 1:num_boot) {
  setTxtProgressBar(prog_bar, value=boot)
  w <- sample.int(n=nrow(draws_mat), size=nrow(draws_mat), replace=TRUE)
  sens_mat_boot <- GetSensitivityFromGrads(grad_mat[, w, drop=FALSE], draws_mat[w, , drop=FALSE])
  sens_mat_array[boot,,] <- sens_mat_boot
  sens_mat_normalized_array[boot,,] <-
    NormalizeSensitivityMatrix(sens_mat_boot, draws_mat[w, , drop=FALSE])
}
close(prog_bar)


GetArrayQuantiles <- function(sens_mat_array, alpha=0.05) {
  lower <- apply(sens_mat_array, MARGIN=c(2, 3),
                 function(x) { quantile(x, alpha) })
  upper <- apply(sens_mat_array, MARGIN=c(2, 3),
                 function(x) { quantile(x, 1 - alpha) })
  
  rownames(lower) <- rownames(upper) <- rownames(sens_mat_list[[1]])
  colnames(lower) <- colnames(upper) <- colnames(sens_mat_list[[1]])
  return(list(lower=lower, upper=upper))
}

sens_mat_quantiles <- GetArrayQuantiles(sens_mat_array)
sens_mat_normalized_quantiles <- GetArrayQuantiles(sens_mat_normalized_array)

sens_mat_quantiles$lower["y_var", ]
sens_mat_quantiles$upper["y_var", ]
sens_mat["y_var", ]

