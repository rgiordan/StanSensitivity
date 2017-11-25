library(mcmcse)
library(coda)
library(doParallel)
registerDoParallel(cores=7)

DrawSamples <- function(num_draws) {
  num_samples <- num_warmup_samples <- num_draws
  sampling_result <- sampling(
    model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples),
    verbose=F, refresh=-1)
  draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
  return(draws_mat)  
}

ResultsFromSamples <- function(par, draws_mat) {
  par_draws <- draws_mat[, par]
  par_mat <- cbind(par_draws, par_draws^2)
  se_results <- mcmcse::mcse.multi(par_mat)
  se_mat <- se_results$cov
  
  # Variance estimate of variance using the Delta method
  # S = 1/n \sum x^2
  # mu = 1/2 \sum x
  # g(S, mu) = S - mu^2
  gbar <- colMeans(par_mat)
  grad_g <- c(-2 * gbar[1], 1)
  g_sd <- sqrt(t(grad_g) %*% se_mat %*% grad_g / length(par_draws))
  
  # Variance estimate with mean held fixed
  g_draws <- (par_draws - mean(par_draws)) ^ 2
  n_eff <-coda::effectiveSize(g_draws)
  g_univariate_sd <- sd(g_draws) / sqrt(n_eff)
  
  return(data.frame(mean=mean(par_draws), var=var(par_draws),
                    var_se=g_sd, var_univariate_se=g_univariate_sd))
}

num_sims <- 200
results_list <- foreach(sim=1:num_sims) %dopar% {
  cat(sim)
  ResultsFromSamples(par="mu", DrawSamples(2000)) %>% mutate(sim=sim)
}

results <- do.call(rbind, results_list)
sd(results$var)
median(results$var_se) / sd(results$var)
median(results$var_univariate_se) / sd(results$var)
sd(results$var_se)
sd(results$var_univariate_se)

hist(results$var)



##################
# Investigate 
draws_mat <- DrawSamples(2000)












###################
# Sanity checking

head(draws_mat)

par <- "mu"

par_draws <- draws_mat[, par]

# Univariate sanity check
par_se <- mcmcse::mcse(par_draws)$se
n_eff <-coda::effectiveSize(par_draws)
n <- length(par_draws)

# About the same I hope
sd(par_draws) / sqrt(n_eff)
par_se

par_se / (sd(par_draws) / sqrt(n))
sqrt(n / n_eff)

# Multivariate

par_mat <- cbind(par_draws, par_draws^2)
cov_mat <- cov(par_mat)
var(par_draws)
se_results <- mcmcse::mcse.multi(par_mat)
se_mat <- se_results$cov

n / n_eff
diag(se_mat) / diag(cov_mat)
