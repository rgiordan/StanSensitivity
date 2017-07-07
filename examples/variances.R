# Compute the variance of the covariances and their ratios.

library(devtools)
# install_local("/home/rgiordan/Documents/git_repos/StanSensitivity/rstansensitivity", force=TRUE)

library(ggplot2)
library(dplyr)
library(reshape2)

boot_results <- BootstrapSensitivityMatrix(draws_mat, sens_result$grad_mat, num_boot=100)

sens_mat_quantiles <-
  GetArrayQuantiles(boot_results$sens_mat_array, alpha=0.1)
sens_mat_normalized_quantiles <-
  GetArrayQuantiles(boot_results$sens_mat_normalized_array, alpha=0.1)

sens_mat_quantiles$lower["y_var", ]
sens_mat_quantiles$upper["y_var", ]
sens_mat["y_var", ]

SensitivityMatrixToDataframe <- function(
    sens_mat, hyperparameter_names, parameter_names) {
  colnames(sens_mat) <- parameter_names
  return(data.frame(sens_mat) %>%
    mutate(hyperparameter=hyperparameter_names) %>%
    melt(id.var="hyperparameter") %>%
    rename(parameter=variable))
}


SummarizeSensitivityMatrices <- function(sens_mat, sens_mat_lower, sens_mat_upper) {
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


sens_df <- SummarizeSensitivityMatrices(
  sens_mat, sens_mat_quantiles$lower, sens_mat_quantiles$upper)
sens_norm_df <- SummarizeSensitivityMatrices(
  sens_mat_normalized,
  sens_mat_normalized_quantiles$lower,
  sens_mat_normalized_quantiles$upper)



ggplot(filter(sens_norm_df, !grepl("weight", hyperparameter))) +
  geom_bar(aes(x=parameter, y=mean_sensitivity, fill=hyperparameter),
           stat="identity", position="dodge") +
  geom_errorbar(aes(x=parameter, ymin=lower_sensitivity,
                    ymax=upper_sensitivity, group=hyperparameter),
                position=position_dodge(0.9), width=0.2)
  
