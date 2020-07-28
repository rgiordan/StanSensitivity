library(testthat)
library(rstansensitivity)
library(rstanarm)
library(sandwich)

context("rstansensitivity")


TestIJ <- function(misspecified, grouped) {
    print(sprintf("Misspecified = %s, grouped = %s", misspecified, grouped))
    set.seed(42)
    num_sims <- 2000
    x1 <- rnorm(num_sims)
    x2 <- rnorm(num_sims)
    if (grouped) {
        z <- sample(floor(num_sims / 10), replace=TRUE, size=num_sims)
    }
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
    if (grouped) {
        loglik_draws <- GroupLogLikelihoodDraws(loglik_draws, z)
    }
    ij_cov <- ComputeIJCovariance(loglik_draws, param_draws)
    scale_mat <- diag(1 / sqrt(diag(ij_cov)))
    rownames(scale_mat) <- colnames(scale_mat) <- rownames(ij_cov)
    
    num_blocks <- 200
    GetBayesIJDifference <- function(loglik_draws, param_draws) {
        ij_cov <- ComputeIJCovariance(loglik_draws, param_draws)
        bayes_cov <- cov(param_draws, param_draws)
        return(ij_cov - ncol(loglik_draws) * bayes_cov)
    }
    ij_bayes_diff_se_list <- GetBlockBootstrapStandardErrors(
        loglik_draws, param_draws, num_blocks=num_blocks, num_draws=50,
        cov_fun=GetBayesIJDifference, 
        show_progress_bar=TRUE)
    
    cc_draws <- ij_bayes_diff_se_list$cov_samples[, 1, 1]
    qqnorm(cc_draws)
    sd(cc_draws)
    ij_bayes_diff_se_list$cov_se[1, 1]
    
    set.seed(42)
    ij_se_list <- GetBlockBootstrapStandardErrors(
        loglik_draws, param_draws, num_blocks=num_blocks, num_draws=50,
        cov_fun=ComputeIJCovariance, 
        show_progress_bar=TRUE)

    if (!misspecified) {
        # Check that IJ matches Bayes if correctly specified.
        # The z score actually rejects quite easily.
        diff <- GetBayesIJDifference(loglik_draws, param_draws)
        z_score <- diff / ij_bayes_diff_se_list$cov_se
        expect_true(max(abs(z_score)) < 4.0, "bayes test")
        
        rel_err <- scale_mat %*% diff %*% scale_mat
        expect_true(max(abs(rel_err)) < 0.2, "Bayes relative error test")
    }

    # Check that IJ matches OLS
    # The z score actually rejects quite easily.
    lm_fit <- lm(y ~ 1 + x1 + x2, df)
    lm_cov <- sandwich::vcovCL(lm_fit, type="HC0", cadjust=FALSE)
    lm_pars <- c("(Intercept)", "x1", "x2")
    
    lm_err <- ncol(loglik_draws) * lm_cov - ij_cov[lm_pars, lm_pars]
    z_score <- lm_err / ij_se_list$cov_se[lm_pars, lm_pars]
    expect_true(max(abs(z_score)) < 4.0, "frequentist z test")
    
    rel_err <- scale_mat[lm_pars, lm_pars] %*% lm_err %*% scale_mat[lm_pars, lm_pars]
    expect_true(max(abs(rel_err)) < 0.2, "frequentist relative error test")
}


test_that("ij_works", {
    #misspecified <- FALSE; grouped <- FALSE
    TestIJ(misspecified=TRUE, grouped=TRUE)
    TestIJ(misspecified=FALSE, grouped=FALSE)
    TestIJ(misspecified=TRUE, grouped=FALSE)
})
