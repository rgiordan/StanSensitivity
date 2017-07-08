library(ggplot2)
library(dplyr)
library(reshape2)


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


BootstrapSensitivityMatrix <- function(
    draws_mat, grad_mat, alpha=0.1, num_boot=200) {
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


SensitivityMatrixToDataframe <- function(
    sens_mat, hyperparameter_names, parameter_names) {
  colnames(sens_mat) <- parameter_names
  return(data.frame(sens_mat) %>%
    mutate(hyperparameter=hyperparameter_names) %>%
    melt(id.var="hyperparameter") %>%
    rename(parameter=variable))
}


SummarizeSensitivityMatrices <- function(
    sens_mat, sens_mat_lower, sens_mat_upper) {
  sens_df <- rbind(
    SensitivityMatrixToDataframe(
      sens_mat,
      hyperparameter_names=rownames(sens_mat),
      parameter_names=colnames(sens_mat)) %>%
      mutate(measure="mean_sensitivity"),
    SensitivityMatrixToDataframe(
      sens_mat_upper,
      hyperparameter_names=rownames(sens_mat),
      parameter_names=colnames(sens_mat)) %>%
      mutate(measure="upper_sensitivity"),
    SensitivityMatrixToDataframe(
      sens_mat_lower,
      hyperparameter_names=rownames(sens_mat),
      parameter_names=colnames(sens_mat)) %>%
      mutate(measure="lower_sensitivity")) %>%
    dcast(hyperparameter + parameter ~ measure) %>%
    filter(parameter != "lp__")
  return(sens_df)
}


GetTidyResult <- function(draws_mat, sens_result, num_boot=500, alpha=0.05) {
  boot_results <- BootstrapSensitivityMatrix(
    draws_mat, sens_result$grad_mat, num_boot=100)
  sens_mat_quantiles <-
    GetArrayQuantiles(boot_results$sens_mat_array, alpha=0.1)
  sens_mat_normalized_quantiles <-
    GetArrayQuantiles(boot_results$sens_mat_normalized_array, alpha=0.1)
  sens_df <- SummarizeSensitivityMatrices(
    sens_mat, sens_mat_quantiles$lower, sens_mat_quantiles$upper)
  sens_norm_df <- SummarizeSensitivityMatrices(
    sens_mat_normalized,
    sens_mat_normalized_quantiles$lower,
    sens_mat_normalized_quantiles$upper)
  return(list(sens_df=sens_df, sens_norm_df=sens_norm_df))
}
