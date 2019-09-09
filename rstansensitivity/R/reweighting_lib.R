# These are experimental ways of doing a different design.  i'm not going to
# export any for the moment, but I'll keep them here as a placeholder.

library(tidybayes)
library(dplyr)


GetWeightMatrixSecondDerivative <- function(stanfit, log_lik_name="log_lik") {
  log_lik <- t(as.matrix(stanfit, log_lik_name))
  log_lik <- log_lik - rowMeans(log_lik)

  SecondDerivative <- function(w1, w2, draws, base_w=1) {
    log_lik1 <- t(w1 - base_w) %*% log_lik
    log_lik2 <- t(w2 - base_w) %*% log_lik
    d2dw_mat <- log_lik1 * log_lik2
    d2dw_mat <- d2dw_mat - mean(d2dw_mat)
    d2dw_mat %*% draws / nrow(draws)
  }

  return(SecondDerivative)
}


GetWeightMatrixDerivative <- function(stanfit, log_lik_name="log_lik") {
    log_lik <- t(as.matrix(stanfit, log_lik_name))
    log_lik <- log_lik - rowMeans(log_lik)

    Derivative <- function(draws) {
        log_lik %*% draws / nrow(draws)
    }

    return(Derivative)
}


GetWeightMatrixPredictor <- function(stanfit, log_lik_name="log_lik") {
    Derivative <- GetWeightMatrixDerivative(
        stanfit=stanfit, log_lik_name=log_lik_name)
    PredictDiff <- function(w, draws, base_w=1) {
        t(w - base_w) %*% Derivative(draws)
    }

    return(PredictDiff)
}


GetTidyWeightPredictor <- function(...,
                                   stanfit=NULL,
                                   PredictDiff=NULL,
                                   draws=NULL,
                                   base_w=1) {
    pars <- enquos(...)
    if (!xor(is.null(stanfit), is.null(PredictDiff))) {
        stop("You must specify `stanfit` or `PredictDiff` but not both.")
    }
    if (!xor(is.null(stanfit), is.null(draws))) {
        stop("You must specify `stanfit` or `draws` but not both.")
    }
    if (!is.null(stanfit)) {
        PredictDiff <- GetWeightMatrixPredictor(stanfit)
        # Take the quoted variables, convert to strings, and strip the index
        # for passing to Stan directly,
        par_names <-
          lapply(pars, as_label) %>%
          sapply(function(x) { sub("\\[.*$", "", x) })
        draws <- as.matrix(stanfit, par_names)
    }
    par_means <-
        gather_draws(t(colMeans(draws)), !!!pars) %>%
        RemoveExtraTidyColumns() %>%
        rename(base_mean=.value)
    join_vars <- setdiff(names(par_means), "base_mean")

    Predict <- function(w) {
        pred_mat <- PredictDiff(w=w, draws=draws, base_w=base_w)
        pred_df <-
          inner_join(
            gather_draws(pred_mat, !!!pars) %>%
              rename(pred_diff=.value) %>%
              RemoveExtraTidyColumns(),
            par_means,
            by=join_vars) %>%
          mutate(pred_mean=base_mean + pred_diff)
        return(pred_df)
    }
    return(Predict)
}
