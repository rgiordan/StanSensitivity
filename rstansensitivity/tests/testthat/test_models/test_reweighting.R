library(rstan)
library(rstansensitivity)
library(testthat)
rstan_options(auto_write=TRUE)

test_that("exact_derivatives_correct") {
    set.seed(42)
    model <- rstan::stan_model("test_models/reweighting_test.stan")

    n_obs <- 30
    x_sd <- 10
    stan_data <- list(
      n_obs=n_obs,
      x=rnorm(n_obs, mean=0, sd=x_sd),
      x_sd=x_sd,
      mu_sd=1,
      mu_mean=0,
      w=rep(1, n_obs)
    )

    num_mcmc_samples <- 80000
    chains <- 1
    stanfit <- sampling(model, data=stan_data,
                        chains=chains, iter=num_mcmc_samples)

    GetMuDerivs <- function(stan_data) {
      # Closed-form derivatives of E[mu | x, w] with respect to w.
      x_info <- 1 / stan_data$x_sd ^ 2
      prior_info <- 1 / stan_data$mu_sd ^ 2
      num <- x_info * sum(stan_data$x * stan_data$w) +
             prior_info * stan_data$mu_mean
      den <- x_info * sum(stan_data$w) + prior_info
      num_deriv <- x_info * stan_data$x
      den_deriv <- x_info * rep(1, stan_data$n_obs)

      muhat <- num / den
      dmu_dw <- num_deriv / den - muhat * den_deriv / den
      d2mu_dw2 <-
        -(outer(num_deriv, den_deriv) +
          outer(den_deriv, num_deriv)) / (den ^ 2) +
        2 * outer(den_deriv, den_deriv) / (den ^ 3)
      return(list(
        muhat=muhat,
        dmu_dw=dmu_dw,
        d2mu_dw2=d2mu_dw2
      ))
    }

    mu_derivs <- GetMuDerivs(stan_data)
    tidy_summary <- GetTidySummary(stanfit, mu, spread=TRUE)
    stopifnot(abs(tidy_summary$mean - mu_derivs$muhat) /
              tidy_summary$se_mean < 2.5)
    Derivative <- GetWeightMatrixDerivative(stanfit)
    mcmc_dmu_dw <- Derivative(as.matrix(stanfit, "mu"))
    expect_equal(mu_derivs$dmu_dw, mcmc_dmu_dw, tol=1e-6)

    # Test some randomly chosen second derivatives.
    SecondDerivative <- GetWeightMatrixSecondDerivative(stanfit)
    LeaveKOutError <- function(k=floor(0.2 * n_obs)) {
      # Note that if sum(w) == n_obs, then there is no second derivative!
      w_new <- rep(1, n_obs)
      w_new[sample(n_obs, k, replace=FALSE)]  <- 0
      d2_mcmc <- SecondDerivative(w_new, w_new, mu_draws)
      d2 <- t(w_new - 1) %*% mu_derivs$d2mu_dw2 %*% (w_new - 1)
      return(data.frame(d2_mcmc=as.numeric(d2_mcmc), d2=d2))
    }

    d2_df <-
      do.call(bind_rows, lapply(1:100, function(x) { LeaveKOutError() })) %>%
      mutate(err=d2_mcmc - d2, rel_error=abs(err) / max(abs(d2)))
    expect_true(max(d2_df$rel_error) < 0.15)
}



test_that("multiple_parameters_work") {
    set.seed(42)

    model <- rstan::stan_model("test_models/reweighting2_test.stan")

    n_obs <- 100
    x_dim <- 3
    mu_true <- runif(x_dim)
    sigma_true <- 0.4
    stan_data <- list(
      n_obs=n_obs,
      x_dim=x_dim,
      x=(matrix(rnorm(n_obs * x_dim), nrow=n_obs) + mu_true) * sigma_true,
      w=rep(1, n_obs)
    )

    num_mcmc_samples <- 3000
    chains <- 1
    stanfit <- sampling(model, data=stan_data,
                        chains=chains, iter=num_mcmc_samples)

    stan_data_w <- stan_data
    stan_data_w$w <- as.numeric(rmultinom(1, n_obs, rep(1 / n_obs, n_obs)))
    stanfit_w <- sampling(model, data=stan_data_w,
                          chains=chains, iter=num_mcmc_samples)
    mcmc_df <-
      gather_draws(stanfit_w, mu[i], sigma) %>%
      group_by(i, .variable) %>%
      summarize(mean=mean(.value), draw_sd=sd(.value))

    TidyPredict <- GetTidyWeightPredictor(stanfit=stanfit, mu[i], log_sigma)

    GetPredDf <- function(TidyPredict) {
      inner_join(
        TidyPredict(stan_data_w$w),
        mcmc_df,
        by=c("i", ".variable")) %>%
        mutate(error=pred_mean - mean,
               rel_error=abs(pred_mean - mean) / draw_sd)
    }

    pred_df1 <-
      GetTidyWeightPredictor(
          stanfit=stanfit, mu[i], sigma) %>%
      GetPredDf()

    draws <- as.matrix(stanfit, pars=c("mu", "sigma"))
    MatrixPredict <- GetWeightMatrixPredictor(stanfit)
    pred_df2 <-
      GetTidyWeightPredictor(
          draws=draws, PredictDiff=MatrixPredict, mu[i], sigma) %>%
      GetPredDf()


    tidy_summary <- GetTidySummary(stanfit, mu[i], log_sigma, spread=TRUE)
    tidy_summary_w <- GetTidySummary(stanfit_w, mu[i], log_sigma, spread=TRUE)
    pred_df1 %>%
      inner_join(
        select(tidy_summary, i, .variable, se_mean),
        by=c(".variable", "i")) %>%
      inner_join(
        select(tidy_summary_w, i, .variable, se_mean),
        by=c(".variable", "i"),
        suffix=c("", "_w"))

    expect_true(all(pred_df1 == pred_df2, na.rm=TRUE))
    expect_true(all(pred_df1$rel_error < 0.07))
}
