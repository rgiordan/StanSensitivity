library(ggplot2)
library(dplyr)
library(reshape2)


####################################
# Perturb and re-fit with importance sampling.

GetPosterior <- function(stan_data) {
  result <- sampling(model, data=stan_data, chains=1,
                     iter=(num_samples + num_warmup_samples),
                     verbose=FALSE)
  return(result)
}

perturb_par <- "cauchy_scale_alpha"
epsilon <- 10.0

# Importance sampling.
se_mean <- summary(sampling_result)$summary[, "se_mean"]
min_epsilon <- 2.0 * min(se_mean / abs(sens_result$sens_mat[perturb_par, ]))
if (epsilon < min_epsilon) {
  warning("The expected change is less than twice the mean standard error for every parameter.")
}

post_list <- list()
imp_list <- list()
num_eff_samples_list <- list()
eps_length <- 5
eps_vec <- seq(0, epsilon, length.out=eps_length)
for (this_eps in eps_vec) {

  stan_data_perturb <- stan_data
  stan_data_perturb[[perturb_par]] <- stan_data_perturb[[perturb_par]] + this_eps
  post_list[[length(post_list) + 1]] <- GetPosterior(stan_data_perturb)
  
  # imp_sens_par_list <- stan_sensitivity_list$sens_par_list
  # imp_sens_par_list[[perturb_par]] <- imp_sens_par_list[[perturb_par]] + this_eps
  # 
  # imp_results <- GetImportanceSamplingFromModelFit(
  #   sampling_result, draws_mat, stan_sensitivity_list,
  #   imp_sens_par_list, lp_vec=sens_result$lp_vec)
  # 
  # imp_diff <- imp_results$imp_means - colMeans(draws_mat)
  # num_eff_samples_list[[length(num_eff_samples_list) + 1]] <-
  #   imp_results$eff_num_imp_samples
  # imp_list[[length(imp_list) + 1]] <-
  #   data.frame(t(imp_diff)) %>%
  #   mutate(epsilon=this_eps) %>%
  #   melt(id.vars="epsilon") %>%
  #   rename(parameter=variable, imp_diff=value)
  
}


MakeTidyMCMCList <- function(post_list, colname) {
  post_mat <- do.call(
    rbind, lapply(post_list, function(result) { summary(result)$summary[, colname] }))
  post_df <- as_tibble(post_mat)
  colnames(post_df) <- colnames(post_mat)
  post_df$epsilon <- eps_vec
  post_df <-
    melt(post_df, id.vars="epsilon") %>%
    rename(parameter=variable)
  post_df[colname] <- post_df$value
  post_df <- select(post_df, -value)
  return(post_df)
}


post_mean_df <- MakeTidyMCMCList(post_list, "mean")
post_se_df <- MakeTidyMCMCList(post_list, "se_mean")
post_df <- inner_join(post_mean_df, post_se_df, by=c("epsilon", "parameter"))

comparison_df <- inner_join(
  post_df, filter(tidy_results$sens_df, hyperparameter==perturb_par), by="parameter")
comparison_df <- inner_join(
  comparison_df, filter(comparison_df, epsilon == 0.0) %>%
                 select("parameter", "mean", "se_mean"),
  by="parameter", suffix=c("", "_original"))
comparison_df <-
  mutate(comparison_df,
         true_sensitivity=mean - mean_original,
         true_sensitivity_se=sqrt(se_mean^2 + se_mean_original^2)) %>%
  filter(epsilon != 0.0)

ggplot(comparison_df) +
  geom_point(aes(y=epsilon * mean_sensitivity, x=true_sensitivity, color=parameter), size=3) +
  geom_errorbarh(aes(y=epsilon * mean_sensitivity,
                     x=true_sensitivity,
                     xmin=true_sensitivity - 2 * true_sensitivity_se,
                     xmax=true_sensitivity + 2 * true_sensitivity_se,
                    color=parameter)) +
  geom_abline(aes(slope=1, intercept=0))




max_eps <- max(unique(imp_df$epsilon)[unlist(num_eff_samples_list) > 500])

ggplot(filter(imp_df, parameter %in% sensitive_params[1:10])) +
  geom_point(aes(x=epsilon, y=imp_diff, color=parameter)) +
  geom_line(aes(x=epsilon, y=mean_sensitivity * epsilon, color=parameter)) +
  geom_vline(aes(xintercept=max_eps))



# Re-run MCMC
stan_data_perturb <- stan_data
stan_data_perturb[[perturb_par]] <- stan_data_perturb[[perturb_par]] + epsilon
sampling_result_perturb <- sampling(model, data=stan_data_perturb, chains=1,
                                    iter=(num_samples + num_warmup_samples))

draws_mat_perturb <- extract(sampling_result_perturb, permute=FALSE)[,1,]
perturb_se <-
  data.frame(t(summary(sampling_result_perturb)$summary[, "se_mean"])) %>%
  mutate(method="mcmc") %>%
  melt(id.vars="method") %>%
  rename(parameter=variable, se=value)

mcmc_diff <- colMeans(draws_mat_perturb) - colMeans(draws_mat)
mcmc_diff_df <-
  data.frame(t(mcmc_diff)) %>%
  mutate(method="mcmc") %>%
  melt(id.vars=c("method")) %>%
  rename(parameter=variable, mean_diff=value) %>%
  inner_join(filter(tidy_results$sens_df, hyperparameter==perturb_par), by="parameter") %>%
  inner_join(perturb_se, by=c("parameter", "method"))

mcmc_diff_df <-
  mutate(mcmc_diff_df,
         top_change=abs(mean_diff) > quantile(abs(mean_diff), 0.9))

ggplot(filter(mcmc_diff_df, parameter %in% sensitive_params)) +
  geom_point(aes(y=mean_diff, x=mean_sensitivity * epsilon, color=parameter), size=3) +
  geom_errorbar(aes(ymin=mean_diff - 2 * se, ymax=mean_diff + 2 * se, x=mean_sensitivity * epsilon)) +
  geom_errorbarh(aes(y=mean_diff, x=mean_sensitivity * epsilon,
                     xmin=lower_sensitivity * epsilon, xmax=upper_sensitivity * epsilon)) +
  geom_abline(aes(slope=1, intercept=0))







############# Old:

# This is intended to be run after run_examples.
sens_mat <- sens_result$sens_mat

GetPosterior <- function(stan_data) {
  result <- sampling(model, data=stan_data, chains=1,
                     iter=(num_samples + num_warmup_samples), verbose=FALSE)
  return(summary(result)$summary)
}

# hyperparam_rows <- !grepl("weights", rownames(sens_mat))
hyperparam_names <- rownames(sens_mat)
param_means <- summary(sampling_result)$summary[, "mean"]

pert_mat <- matrix(NA, length(hyperparam_names), ncol(sens_mat))
rownames(pert_mat) <- hyperparam_names
colnames(pert_mat) <- colnames(sens_mat)
param_names <- setdiff(colnames(pert_mat), "lp__")

# Set the perturbation well outside the smallest standard error range.
epsilon <- 2 * max(summary(sampling_result)$summary[param_names, "se_mean"])

for (i in 1:length(hyperparam_names)) {
  hp_name <- hyperparam_names[i]
  print(sprintf('Perturbing %s', hp_name))
  stan_data_perturb <- stan_data
  stan_data_perturb[[hp_name]] <- stan_data_perturb[[hp_name]]  + epsilon
  pert_mat[i, ] <- GetPosterior(stan_data_perturb)[, "mean"] - param_means
}

pert_mat[, param_names]
sens_mat[, param_names] * epsilon



