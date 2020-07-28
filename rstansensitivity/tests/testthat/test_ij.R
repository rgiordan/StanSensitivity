library(testthat)
library(rstansensitivity)
library(rstanarm)
library(sandwich)

context("rstansensitivity")


TestIJ <- function(misspecified) {
    set.seed(42)
    num_sims <- 5000
    x1 <- rnorm(num_sims)
    x2 <- rnorm(num_sims)
    sigma2 <- 3 * (x1^2 + x2^2)
    
    if (misspecified) {
        eps <- rnorm(num_sims, sd=sqrt(sigma2))
    } else {
        eps <- rnorm(num_sims, sd=median(sqrt(sigma2)))
    }
    #plot(sqrt(x1^2 + x2^2), eps)
    df <- data.frame(y=eps, x1=x1, x2=x2)
    
    rstan_fit <- rstanarm::stan_glm(y ~ 1 + x1 + x2, df, family=gaussian())
    param_draws <- as.matrix(rstan_fit)
    param_draws[, "sigma"] <- log(param_draws[, "sigma"])
    bayes_cov <- cov(param_draws, param_draws)
    
    loglik_draws <- log_lik(rstan_fit)
    ij_cov <- ComputeIJCovariance(loglik_draws, param_draws)
    
    GetBayesIJDifference <- function(loglik_draws, param_draws) {
        ij_cov <- ComputeIJCovariance(loglik_draws, param_draws)
        bayes_cov <- cov(param_draws, param_draws)
        return(ij_cov - ncol(loglik_draws) * bayes_cov)
    }
    ij_bayes_diff_se_list <- GetBlockBootstrapStandardErrors(
        loglik_draws, param_draws, num_blocks=100, num_draws=20,
        cov_fun=GetBayesIJDifference, 
        show_progress_bar=TRUE)
    
    ij_se_list <- GetBlockBootstrapStandardErrors(
        loglik_draws, param_draws, num_blocks=100, num_draws=20,
        cov_fun=ComputeIJCovariance, 
        show_progress_bar=TRUE)
    
    if (!misspecified) {
        # Check that IJ matches Bayes if correctly specified.
        # This is just a coarse check that does not account for multiplicity.
        diff <- GetBayesIJDifference(loglik_draws, param_draws)
        z_score <- diff / ij_bayes_diff_se_list$cov_se
        expect_true(max(abs(z_score)) < 2.5)
    }

    # Check that IJ matches OLS
    # This is just a coarse check that does not account for multiplicity.
    lm_fit <- lm(y ~ 1 + x1 + x2, df)
    lm_cov <- sandwich::vcovHC(lm_fit)
    lm_pars <- c("(Intercept)", "x1", "x2")
    
    num_sims * lm_cov
    ij_cov[lm_pars, lm_pars]
    lm_err <- num_sims * lm_cov - ij_cov[lm_pars, lm_pars]
    z_score <- lm_err / ij_se_list$cov_se[lm_pars, lm_pars]
    expect_true(max(abs(z_score)) < 2.5)
}


test_that("ij_works", {
    TestIJ(TRUE)
    TestIJ(FALSE)
})
