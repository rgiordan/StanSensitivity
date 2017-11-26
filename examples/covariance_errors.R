library(mcmcse)
library(coda)
library(doParallel)
registerDoParallel(cores=7)

DrawSamples <- function(num_draws, draws_only=TRUE) {
  num_samples <- num_warmup_samples <- num_draws
  sampling_result <- sampling(
    model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples),
    verbose=F, refresh=-1)
  if (draws_only) {
    draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
    return(draws_mat)  
  } else {
    return(sampling_result)
  }
}


#####################
# Standard errors of covariance

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


#####################
# Standard errors of correlation (not useful)

GetCorrelationGradient <- function(x_draws, y_draws) {
  x_mean <- mean(x_draws)
  y_mean <- mean(y_draws)
  sigma_xy <- mean(x_draws * y_draws) - x_mean * y_mean
  sigma_x <- sqrt(mean(x_draws^2) - x_mean^2)
  sigma_y <- sqrt(mean(y_draws^2) - y_mean^2)
  
  arg_draws <- cbind(x_draws * y_draws, x_draws^2, y_draws^2, x_draws, y_draws)
  arg_cov_mat <- mcmcse::mcse.multi(arg_draws)$cov
  grad_g <- c(1 / (sigma_x * sigma_y),
              -0.5 * sigma_xy / (sigma_x^3 * sigma_y),
              -0.5 * sigma_xy / (sigma_x * sigma_y^3),
              (sigma_xy * x_mean / (sigma_x^2) - y_mean) / (sigma_x * sigma_y),
              (sigma_xy * y_mean / (sigma_y^2) - x_mean) / (sigma_x * sigma_y))
  return(grad_g)
}

GetCorrelationSE <- function(x_draws, y_draws) {
  arg_draws <- cbind(x_draws * y_draws, x_draws^2, y_draws^2, x_draws, y_draws)
  arg_cov_mat <- mcmcse::mcse.multi(arg_draws)$cov
  grad_g <- GetCorrelationGradient(x_draws, y_draws)
  g_se <- as.numeric(sqrt(t(grad_g) %*% arg_cov_mat %*% grad_g / nrow(arg_draws)))

  return(g_se)
}

# Numerically check the gradient.
library(numDeriv)

PackCorrelationPar <- function(x_draws, y_draws) {
  return(c(
    mean(x_draws * y_draws),
    mean(x_draws^2),
    mean(y_draws^2),
    mean(x_draws),
    mean(y_draws)
  ))
}

GetCorrelation <- function(par) {
  xy_mean <- par[1]
  x2_mean <- par[2]
  y2_mean <- par[3]
  x_mean <- par[4]
  y_mean <- par[5]
  
  sigma_x <- sqrt(x2_mean - x_mean^2)
  sigma_y <- sqrt(y2_mean - y_mean^2)
  sigma_xy <- xy_mean - x_mean * y_mean
  return(sigma_xy / (sigma_x * sigma_y))
}

x_draws <- grad_mat[1, ]
y_draws <- draws_mat[, 1]
g_grad_numeric <- grad(GetCorrelation,
                       PackCorrelation(x_draws, y_draws))
g_grad <- GetCorrelationGradient(x_draws, y_draws)
g_grad / g_grad_numeric


arg_draws <- cbind(x_draws * y_draws, x_draws^2, y_draws^2, x_draws, y_draws)
arg_sd_vec <- apply(arg_draws, MARGIN=2, sd)
arg_sd_mat <- diag(arg_sd_vec)
arg_draws_norm <- arg_draws / rep(arg_sd_vec, each=nrow(arg_draws))


arg_cov_mat <- mcmcse::mcse.multi(arg_draws)$cov
arg_cov_mat_norm <- mcmcse::mcse.multi(arg_draws_norm)$cov

# It's scale invariant
(arg_sd_mat %*% arg_cov_mat_norm %*% arg_sd_mat) / arg_cov_mat

grad_g <- GetNormalizedCovarianceGradient(x_draws, y_draws)
g_se <- as.numeric(sqrt(t(grad_g) %*% arg_cov_mat %*% grad_g / nrow(arg_draws)))


#####################
# Standard errors of normalized covariance

PackNormalizedCovariancePar <- function(x_draws, y_draws) {
  return(cbind(
    x_draws * y_draws,
    x_draws ^ 2,
    x_draws,
    y_draws
  ))
}

GetNormalizedCovarianceGradient <- function(par) {
  mean_xy <- par[1]
  mean_x2 <- par[2]
  mean_x <- par[3]
  mean_y <- par[4]
  sigma_xy <- mean_xy - mean_x * mean_y
  sigma_x <- sqrt(mean_x2 - mean_x^2)
  return(c(
    1 / sigma_x,
    -0.5 * sigma_xy / (sigma_x^3),
    (sigma_xy * mean_x / (sigma_x^2) - mean_y) / sigma_x,
    -mean_x / sigma_x
  ))
}

GetNormalizedCovariance <- function(par) {
  mean_xy <- par[1]
  mean_x2 <- par[2]
  mean_x <- par[3]
  mean_y <- par[4]
  sigma_xy <- mean_xy - mean_x * mean_y
  sigma_x <- sqrt(mean_x2 - mean_x^2)
  return(sigma_xy / sigma_x)  
}

GetNormalizedCovarianceSE <- function(x_draws, y_draws) {
  par_draws <- PackNormalizedCovariancePar(x_draws, y_draws)
  par_cov_mat <- mcmcse::mcse.multi(par_draws)$cov
  grad_g <- GetNormalizedCovarianceGradient(colMeans(par_draws))
  return(
    as.numeric(sqrt(t(grad_g) %*% par_cov_mat %*% grad_g / nrow(par_draws))))
}

x_draws <- grad_mat[1, ]
y_draws <- draws_mat[, 1]
par <- colMeans(PackNormalizedCovariancePar(x_draws, y_draws))
g_grad_numeric <- grad(GetNormalizedCovariance, par)
g_grad <- GetNormalizedCovarianceGradient(par)
g_grad / g_grad_numeric


GetNormalizedCovarianceSE(x_draws, y_draws)

###################
# Compare with simulations

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



####################
# Sensitivity covariances

grad_mat <- sens_result$grad_mat
head(t(grad_mat))

sens_param_names <- rownames(grad_mat)
param_names <- colnames(draws_mat)

GetSensitivityStandardErrors <- function(draws_mat, grad_mat, fix_mean=FALSE, normalized=FALSE) {
  if (normalized && fix_mean) {
    stop("You cannot specify both fix_mean and normalized.")
  }
  cov_se_mat <- matrix(NA, nrow(grad_mat), ncol(draws_mat))
  for (par_ind in 1:length(param_names)) {
    for (grad_ind in 1:length(sens_param_names)) {
      par_draws <- draws_mat[, param_names[par_ind]]
      grad_draws <- grad_mat[sens_param_names[grad_ind], ]
      if (normalized) {
        cov_se_mat[grad_ind, par_ind] <-
          GetNormalizedCovarianceSE(par_draws, grad_draws)
      } else {
        cov_se_mat[grad_ind, par_ind] <-
          GetCovarianceSE(par_draws, grad_draws, fix_mean=fix_mean)
      }
    }
  }
  return(cov_se_mat)
}

num_sims <- 500
sampling_list <- foreach(sim=1:num_sims) %dopar% {
  DrawSamples(2000, draws_only=FALSE)
}

sens_list <- foreach(sim=1:num_sims) %dopar% {
  sampling_result <- sampling_list[[sim]]
  draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
  sens_result <- GetStanSensitivityFromModelFit(
    sampling_result, draws_mat, stan_sensitivity_model)
  return(sens_result)
}

cov_se_list <- foreach(sim=1:num_sims) %dopar% {
  sampling_result <- sampling_list[[sim]]
  sens_result <- sens_list[[sim]]
  draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
  grad_mat <- sens_result$grad_mat 
  return(GetSensitivityStandardErrors(draws_mat, grad_mat, fix_mean=FALSE))
}

cov_se_fix_list <- foreach(sim=1:num_sims) %dopar% {
  sampling_result <- sampling_list[[sim]]
  sens_result <- sens_list[[sim]]
  draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
  grad_mat <- sens_result$grad_mat 
  return(GetSensitivityStandardErrors(draws_mat, grad_mat, fix_mean=TRUE))
}


cov_norm_se_list <- foreach(sim=1:num_sims) %dopar% {
  sampling_result <- sampling_list[[sim]]
  sens_result <- sens_list[[sim]]
  draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
  grad_mat <- sens_result$grad_mat 
  return(GetSensitivityStandardErrors(draws_mat, grad_mat, normalized=TRUE))
}


cov_array <- simplify2array(lapply(sens_list, function(sens_result) { sens_result$sens_mat }))
cov_sd <- apply(cov_array, 1:2, sd)
cov_sd / apply(simplify2array(cov_se_list), 1:2, median)
cov_sd / apply(simplify2array(cov_se_fix_list), 1:2, median)
apply(simplify2array(cov_se_list), 1:2, sd) / apply(simplify2array(cov_se_list), 1:2, median)


cov_norm_array <- simplify2array(lapply(sens_list, function(sens_result) { sens_result$sens_mat_normalized }))
cov_norm_sd <- apply(cov_norm_array, 1:2, sd)
cov_norm_sd / apply(simplify2array(cov_norm_se_list), 1:2, median)
apply(simplify2array(cov_norm_se_list), 1:2, sd) /
  apply(simplify2array(cov_norm_se_list), 1:2, median)


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
