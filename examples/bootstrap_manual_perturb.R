
# This is intended to be run after run_examples and manual_perturb on a model that
# measures sensitivity to data weights.

GetPosterior <- function(stan_data) {
  sampling_result <- sampling(
    model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples))
  return(summary(sampling_result)$summary)  
}

###########################
# Compare with bootstrap sampling.
sens_mat <- sens_result$sens_mat
param_names <- colnames(sens_mat)
param_means <- colMeans(draws_mat)

weight_rows <- grepl("weights", rownames(sens_mat))
weight_sens <- sens_mat[weight_rows, param_names, drop=FALSE]

num_boot <- 200
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

mean_boot_df <-
  rbind(
    data.frame(mean_boot) %>%
      mutate(draw=1:n()) %>%
      melt(id.vars="draw") %>%
      mutate(method="mcmc")
,
  data.frame(mean_boot_lr) %>%
    mutate(draw=1:n()) %>%
    melt(id.vars="draw") %>%
    mutate(method="linear")
  ) %>% filter(variable != "lp__")

boot_df <-
  dcast(mean_boot_df, draw + variable ~ method) %>%
  filter(variable != "lp__")

ggplot(boot_df) +
  geom_point(aes(x=mcmc, y=linear)) +
  geom_abline(aes(slope=1, intercept=0), color="gray") +
  facet_grid(~ variable) +
  xlab("Bootstrapped posterior mean") +
  ylab("Linear approximate bootstrapped posterior mean")

ggplot(boot_df) +
  geom_density(aes(x=mcmc, color="Full bootstrap"), n=50, lwd=2) +
  geom_density(aes(x=linear, color="Linear bootstrap"), n=50, lwd=2)

method_labeller <- function(method) {
  return(ifelse(method == "linear", "Linear bootstrap", "Full bootstrap"))
}

ggplot(mean_boot_df) +
  geom_histogram(aes(x=value, color=variable), lwd=1) +
  facet_grid(~ method, labeller=as_labeller(method_labeller))

cat("Linear response time: ", as.numeric(lr_time, units="secs"))
cat("Bootstrap time: ", as.numeric(mcmc_time, units="secs"))


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
