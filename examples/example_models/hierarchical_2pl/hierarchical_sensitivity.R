library(rstan)
rstan_options(auto_write=TRUE)
library(ggplot2)
library(dplyr)
library(reshape2)
library(rstansensitivity)


# Set this to be the appropriate location of the repo on your computer.
example_directory <- file.path(
  Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/examples/example_models/hierarchical_2pl")


library(mirt)
sat <- key2binary(SAT12,
                  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))

stan_data <- list(I = ncol(sat),
                 J = nrow(sat),
                 N = length(sat),
                 ii = rep(1:ncol(sat), each = nrow(sat)),
                 jj = rep(1:nrow(sat), times = ncol(sat)),
                 y = as.vector(sat))

# Set the hyperparameters.
stan_data$lkj_concentration <- 4.0
stan_data$mu_loc <- 0.
stan_data$mu_1_scale <-1.
stan_data$mu_2_scale <- 5.
stan_data$tau_loc <- 0.1
stan_data$theta_loc <- 0.
stan_data$theta_scale <- 1.

base_model_name <- file.path(example_directory, "hierarchical_2pl.stan")
python_script <- file.path(Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/python/generate_models.py")
model_name <- GenerateSensitivityFromModel(base_model_name, python_script=python_script)
model <- stan_model(paste(model_name, "_generated.stan", sep=""))

num_warmup_samples <- num_samples <- 1000
mcmc_time <- Sys.time()
sampling_result <- sampling(model, data=stan_data, chains=1, iter=num_warmup_samples + num_samples)
mcmc_time <- Sys.time() - mcmc_time

###################
# Sensitivity

draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
stan_sensitivity_list <- GetStanSensitivityModel(sampling_result, model_name, stan_data)

sens_time <- Sys.time()
sens_result <- GetStanSensitivityFromModelFit(sampling_result, draws_mat, stan_sensitivity_list)
sens_time <- Sys.time()- sens_time


##################################
# Inspect the results.

# Warning: the uncertainty estimates on the sensitivity are currently underestimated,
# as they do not take into account autocorrelation in the MCMC chain.
#library(devtools)
#install_local("/home/rgiordan/Documents/git_repos/StanSensitivity/rstansensitivity/", force=TRUE)

tidy_results <- GetTidyResult(draws_mat, sens_result)
sens_norm_df <- tidy_results$sens_norm_df
sens_norm_df$parameter_base <- sub("\\..*$", "", sens_norm_df$parameter)

primary_parameters <- c("alpha", "beta", "mu", "tau")

ggplot(filter(sens_norm_df, abs(mean_sensitivity) > 1, parameter_base %in% primary_parameters)) +
  geom_bar(aes(x=parameter, y=mean_sensitivity, fill=hyperparameter),
           stat="identity", position="dodge") +
  geom_errorbar(aes(x=parameter, ymin=lower_sensitivity,
                    ymax=upper_sensitivity, group=hyperparameter),
                position=position_dodge(0.9), width=0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Normalized local sensitivity") +
  ggtitle(sprintf("Hyperparameter sensitivity for model %s", sub(".*/", "", model_name)))

sensitive_params <- 
  filter(sens_norm_df, abs(mean_sensitivity) > 1, parameter_base %in% primary_parameters)$parameter

####################################
# Perturb and re-fit


# Importance sampling:

perturb_par <- "theta_scale"
epsilon <- 0.5

# Importance sampling.
se_mean <- summary(sampling_result)$summary[, "se_mean"]
min_epsilon <- 2.0 * min(se_mean / abs(sens_result$sens_mat[perturb_par, ]))
if (epsilon < min_epsilon) {
  warning("The expected change is less than twice the mean standard error for every parameter.")
}

imp_list <- list()
num_eff_samples_list <- list()
eps_length <- 5
for (this_eps in seq(0, epsilon, length.out=eps_length)) {
  imp_sens_par_list <- stan_sensitivity_list$sens_par_list
  imp_sens_par_list[[perturb_par]] <- imp_sens_par_list[[perturb_par]] + this_eps
  
  imp_results <- GetImportanceSamplingFromModelFit(
    sampling_result, draws_mat, stan_sensitivity_list,
    imp_sens_par_list, lp_vec=sens_result$lp_vec)
  
  imp_diff <- imp_results$imp_means - colMeans(draws_mat)
  num_eff_samples_list[[length(num_eff_samples_list) + 1]] <-
    imp_results$eff_num_imp_samples
  imp_list[[length(imp_list) + 1]] <-
    data.frame(t(imp_diff)) %>%
    mutate(epsilon=this_eps) %>%
    melt(id.vars="epsilon") %>%
    rename(parameter=variable, imp_diff=value)
}

imp_df <- do.call(rbind, imp_list) %>%
  inner_join(filter(tidy_results$sens_df, hyperparameter==perturb_par), by="parameter")

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


