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


SensitivityMatrixToDataframe <- function(
    sens_mat, hyperparameter_names, parameter_names) {
  colnames(sens_mat) <- parameter_names
  return(data.frame(sens_mat) %>%
    mutate(hyperparameter=hyperparameter_names) %>%
    melt(id.var="hyperparameter") %>%
    rename(parameter=variable))
}


SummarizeSensitivityMatrices <- function(sens_mat, sens_se) {
    sens_df <- rbind(
      SensitivityMatrixToDataframe(
        sens_mat,
        hyperparameter_names=rownames(sens_mat),
        parameter_names=colnames(sens_mat)) %>%
        mutate(measure="sensitivity_mean"),
      SensitivityMatrixToDataframe(
        sens_se,
        hyperparameter_names=rownames(sens_se),
        parameter_names=colnames(sens_se)) %>%
        mutate(measure="sensitivity_se")) %>%
      dcast(hyperparameter + parameter ~ measure) %>%
      filter(parameter != "lp__")
    return(sens_df)
}


GetTidyResult <- function(draws_mat, sens_result) {
  sens_se <- GetSensitivityStandardErrors(
      draws_mat, sens_result$grad_mat, fix_mean=FALSE, normalized=FALSE)
  norm_sens_se <- GetSensitivityStandardErrors(
      draws_mat, sens_result$grad_mat, fix_mean=FALSE, normalized=TRUE)
  sens_df <- SummarizeSensitivityMatrices(sens_result$sens_mat, sens_se)
  sens_norm_df <- SummarizeSensitivityMatrices(
    sens_result$sens_mat_normalized, norm_sens_se)
  return(list(sens_df=sens_df, sens_norm_df=sens_norm_df))
}


PlotSensitivities <- function(sens_df, se_num=2) {
  return(
    ggplot(sens_df) +
      geom_bar(aes(x=parameter, y=sensitivity_mean, fill=hyperparameter),
               stat="identity", position="dodge") +
      geom_errorbar(aes(x=parameter,
                        ymin=sensitivity_mean - se_num  * sensitivity_se,
                        ymax=sensitivity_mean + se_num  * sensitivity_se,
                        group=hyperparameter),
                    position=position_dodge(0.9), width=0.2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_discrete(name="Hyperparameter") +
      ylab("Sensitivity") + xlab("Parameter")
  )
}
