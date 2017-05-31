library(ggplot2)
library(dplyr)
library(reshape2)


# This is intended to be run after run_examples.

GetPosterior <- function(stan_data) {
  result <- sampling(model, data=stan_data, chains=1, iter=num_samples * 2,
                     refresh=-1, verbose=FALSE, show_messages=FALSE)
  return(summary(result)$summary)
}

hyperparam_rows <- !grepl("weights", rownames(sens_mat))
hyperparam_names <- rownames(sens_mat)[param_rows]
param_means <- summary(result)$summary[, "mean"]

pert_mat <- matrix(NA, length(hyperparam_names), ncol(sens_mat))
rownames(pert_mat) <- hyperparam_names
colnames(pert_mat) <- colnames(sens_mat)
param_names <- setdiff(colnames(pert_mat), "lp__")

# Set the perturbation well outside the standard error range.
epsilon <- 10 * min(summary(result)$summary[param_names, "se_mean"])

for (i in 1:length(hyperparam_names)) {
  hp_name <- hyperparam_names[i]
  stan_data_perturb <- stan_data
  stan_data_perturb[[hp_name]] <- stan_data_perturb[[hp_name]]  + epsilon
  pert_mat[i, ] <- GetPosterior(stan_data_perturb)[, "mean"] - param_means
}

pert_mat[, param_names]
sens_mat[hyperparam_names, param_names] * epsilon


###########################
# Compare with bootstrap sampling.

weight_rows <- grepl("weights", rownames(sens_mat))
weight_sens <- sens_mat[weight_rows, param_names, drop=FALSE]

num_boot <- 50
num_samples <- 1000 # Number of MCMC samples

mean_boot_lr <- matrix(NA, num_boot, length(param_names))
mean_boot <- matrix(NA, num_boot, length(param_names))
se_boot <- matrix(NA, num_boot, length(param_names))
colnames(mean_boot_lr) <- colnames(mean_boot) <- colnames(se_boot) <- param_names

lr_time <- 0.0
mcmc_time <- 0.0

for (b in 1:num_boot) {
  # Draw weights
  w <- as.numeric(rmultinom(n=1, size=stan_data$N_observed, prob=rep(1, stan_data$N_observed)))
  stopifnot(sum(w) == stan_data$N_observed)
  stopifnot(mean(w) == 1)

  # "Linear response" bootstrap
  tic <- Sys.time()
  mean_boot_lr[b, ] <- param_means[param_names] + colSums(weight_sens * (w - 1))
  lr_time <- lr_time + Sys.time() - tic

  # MCMC bootstrap
  stan_data_boot <- stan_data
  stan_data_boot$weights <- w
  tic <- Sys.time()
  post <- GetPosterior(stan_data_boot)
  mcmc_time <- mcmc_time + Sys.time() - tic
  mean_boot[b, ] <- post[param_names, "mean"]
  se_boot[b, ] <- post[param_names, "se_mean"]
}

plot(mean_boot, mean_boot_lr); abline(0, 1)
print(mcmc_time)
print(lr_time)


#################
# Compare the linear bootstrap with the posterior

num_boot <- 10000

mean_boot_lr <- matrix(NA, num_boot, length(param_names))
colnames(mean_boot_lr) <- param_names

for (b in 1:num_boot) {
  # Draw weights
  w <- as.numeric(rmultinom(n=1, size=stan_data$N_observed, prob=rep(1, stan_data$N_observed)))
  stopifnot(sum(w) == stan_data$N_observed)
  stopifnot(mean(w) == 1)

  # "Linear response" bootstrap
  mean_boot_lr[b, ] <- param_means[param_names] + colSums(weight_sens * (w - 1))
}

mean_boot_df <-
  data.frame(mean_boot_lr) %>%
  mutate(draw=1:n()) %>%
  melt(id.vars="draw") %>%
  mutate(method="bootstrap")

draws_df <-
  data.frame(draws_mat) %>%
  mutate(draw=1:n()) %>%
  melt(id.vars="draw") %>%
  mutate(method="mcmc") %>%
  filter(variable != "lp__")

boot_df <- rbind(mean_boot_df, draws_df)

ggplot(boot_df) +
  geom_density(aes(x=value, color=method, group=method)) +
  facet_grid( ~ variable)
