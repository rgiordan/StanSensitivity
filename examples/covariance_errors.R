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

GetCovarianceSE <- function(x_draws, y_draws, fix_mean=FALSE) {
  x_mean <- mean(x_draws)
  y_mean <- mean(y_draws)
  if (fix_mean) {
    g_draws <- x_draws * y_draws - x_mean * y_mean
    n_eff <-coda::effectiveSize(g_draws)
    g_se <- sd(g_draws) / sqrt(n_eff)
  } else {
    arg_draws <- cbind(x_draws * y_draws, x_draws, y_draws)
    arg_cov_mat <- mcmcse::mcse.multi(arg_draws)$cov
    grad_g <- c(1, -1 * y_mean, -1 * x_mean)
    g_se <- as.numeric(sqrt(t(grad_g) %*% arg_cov_mat %*% grad_g / nrow(arg_draws)))
  }
  return(g_se)
}

GetWorstCaseCovarianceSE <- function(x_draws, y_draws, fix_mean=FALSE) {
  x_mean <- mean(x_draws)
  y_mean <- mean(y_draws)
  g1_draws <- x_draws * y_draws - x_mean * y_mean
  ######### Stopped here.
  n_eff <-coda::effectiveSize(g_draws)
  g_se <- sd(g_draws) / sqrt(n_eff)
  return(g_se)
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

num_sims <- 500
draws_list <- foreach(sim=1:num_sims) %dopar% {
  DrawSamples(2000)
}

##############
# Variances

results_list <-
  lapply(1:num_sims,
   function(ind) {
     ResultsFromSamples(par="log_beta", draws_list[ind][[1]])
    })

results <- do.call(rbind, results_list)
sd(results$var)
median(results$var_se) / sd(results$var)
median(results$var_univariate_se) / sd(results$var)
sd(results$var_se)
sd(results$var_univariate_se)

hist(results$var)

##############
# Covariances

par1 <- "log_alpha"
par2 <- "log_beta"
GetCovarianceResults <- function(draws_mat, par1, par2) {
  cov <- cov(draws_mat[, par1], draws_mat[, par2])
  cov_se <- GetCovarianceSE(draws_mat[, par1], draws_mat[, par2])
  cov_fixed_se <- GetCovarianceSE(draws_mat[, par1], draws_mat[, par2], fix_mean=TRUE)
  return(data.frame(cov=cov, cov_se=cov_se, cov_fixed_se=cov_fixed_se))
}

cov_results <-
  do.call(rbind,
          lapply(1:num_sims,
                function(ind) { GetCovarianceResults(draws_list[ind][[1]], par1, par2) }))
sd(cov_results$cov)
median(cov_results$cov_se) / sd(cov_results$cov)
median(cov_results$cov_fixed_se) / sd(cov_results$cov)
sd(cov_results$cov_se) / median(cov_results$cov)
sd(cov_results$cov_fixed_se) / median(cov_results$cov)

qqnorm(cov_results$cov)



##################
# Investigate 

draws_mat <- DrawSamples(2000)

par_draws <- draws_mat[, par]
par_mat <- cbind(par_draws, par_draws^2)
se_results <- mcmcse::mcse.multi(par_mat)
se_mat <- se_results$cov
inv_scale_mat <- solve(diag(sqrt(diag(se_mat)))) 
corr_mat <- inv_scale_mat %*% se_mat %*% inv_scale_mat

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
