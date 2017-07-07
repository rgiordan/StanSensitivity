

# Evaluate at draws using the hyperparameters in sens_par_list.
EvaluateAtDraws <- function(sampling_result, draws_mat, stan_sensitivity_list, sens_par_list,
                            compute_grads=FALSE) {

  num_warmup_samples <- attr(sampling_result, "sim")$warmup
  model_sens_fit <- stan_sensitivity_list$model_sens_fit
  sens_param_names <- stan_sensitivity_list$sens_param_names

  # Get the model gradients with respect to the hyperparameters (and parameters).
  num_samples <- nrow(draws_mat)
  lp_vec <- rep(NA, num_samples)
  if (compute_grads) {
    grad_mat <- matrix(NA, num_samples, length(sens_param_names))
  } else {
    grad_mat <- matrix()
  }

  cat("Evaluating model at the MCMC draws.\n")
  prog_bar <- txtProgressBar(min=1, max=num_samples, style=3)
  for (n in 1:num_samples) {
    setTxtProgressBar(prog_bar, value=n)
    par_list <- get_inits(sampling_result, iter=n + num_warmup_samples)[[1]]
    for (par in ls(par_list)) {
      # Note that get_inits is currently broken:
      # https://github.com/stan-dev/rstan/issues/417
      # ...but this has seemed to fix it (so far):
      if (length(dim(sens_par_list[[par]])) >= 2) {
        sens_par_list[[par]] <-
          array(unlist(par_list[[par]]), dim(par_list[[par]]))
      } else {
        sens_par_list[[par]] <- as.numeric(par_list[[par]])
      }
    }
    pars_free <- unconstrain_pars(model_sens_fit, sens_par_list)
    
    if (compute_grads) {
      glp <- grad_log_prob(model_sens_fit, pars_free)
      grad_mat[n, ] <- glp
      lp_vec[n] <- attr(glp, "log_prob")
    } else {
      lp_vec[n] <- log_prob(model_sens_fit, pars_free)
    }
  }
  close(prog_bar)
  return(list(lp_vec=lp_vec, grad_mat=grad_mat))
}

lp_vec <- EvaluateAtDraws(sampling_result, draws_mat, stan_sensitivity_list,
                          stan_sensitivity_list$sens_par_list,
                          compute_grads=FALSE)$lp_vec
epsilon <- 0.01
imp_sens_par_list <- stan_sensitivity_list$sens_par_list
imp_sens_par_list$y_var <- imp_sens_par_list$y_var + epsilon
imp_lp_vec <- EvaluateAtDraws(sampling_result, draws_mat, stan_sensitivity_list,
                              imp_sens_par_list, compute_grads=FALSE)$lp_vec
imp_weights <- exp(imp_lp_vec - lp_vec)
imp_weights <- imp_weights / sum(imp_weights)
hist(log(imp_weights), 100)
abline(v=-log(nrow(draws_mat)))

imp_mu_diff <- colSums(imp_weights * draws_mat[, "mu", drop=FALSE]) -
               colMeans(draws_mat[, "mu", drop=FALSE])
print(imp_mu_diff)
sens_mat["y_var", "mu"] * epsilon


draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
stan_sensitivity_list <- GetStanSensitivityModel(sampling_result, model_name, stan_data)
sens_result <- GetStanSensitivityFromModelFit(
  sampling_result, draws_mat, stan_sensitivity_list, num_warmup_samples=num_warmup_samples)




#############
# Debug

model_sens_fit <- stan_sensitivity_list$model_sens_fit

n <- 1
sens_par_list <- stan_sensitivity_list$sens_par_list
par_list <- get_inits(sampling_result, iter=n + num_warmup_samples)[[1]]

for (par in ls(par_list)) {
  if (length(dim(sens_par_list[[par]])) >= 2) {
    sens_par_list[[par]] <-
      array(unlist(par_list[[par]]), dim(par_list[[par]]))
  } else {
    sens_par_list[[par]] <- as.numeric(par_list[[par]])
  }
}
pars_free <- unconstrain_pars(model_sens_fit, sens_par_list)

# Should be 0
par_list$mu - draws_mat[n, "mu"]
log_prob(stan_sensitivity_list$model_sens_fit, pars_free) - draws_mat[n, "lp__"]
lp_vec[n] - draws_mat[n, "lp__"]
imp_lp_vec[n] - draws_mat[n, "lp__"]





#########################################
# Perturb and re-draw.  It looks like importance sampling is wrong.

stan_data_perturb <- stan_data
stan_data_perturb$y_var <- stan_data_perturb$y_var + epsilon
sampling_result_perturb <- sampling(model, data=stan_data_perturb, chains=1,
                                    iter=(num_samples + num_warmup_samples))
draws_mat_perturb <- extract(sampling_result_perturb, permute=FALSE)[,1,]
mean(draws_mat_perturb[, "mu"]) - mean(draws_mat[, "mu"])
print(imp_mu_diff)
sens_mat["y_var", "mu"] * epsilon


