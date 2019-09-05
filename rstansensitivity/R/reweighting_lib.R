library(tidybayes)
library(dplyr)


GetWeightMatrixPredictor <- function(stanfit, log_lik_name="log_lik") {
    log_lik <- t(as.matrix(stanfit, log_lik_name))
    log_lik <- log_lik - rowMeans(log_lik)

    PredictDiff <- function(w, draws, base_w=1) {
        t(w - base_w) %*% log_lik %*% draws / nrow(draws)
    }

    return(PredictDiff)
}


CleanGather <- function(mat, ...) {
  pars <- enquos(...)
  gather_draws(mat, !!!pars) %>%
    select(-.chain, -.iteration, -.draw)
}


GetTidyWeightPredictor <- function(stanfit=NULL,
                                   PredictDiff=NULL,
                                   draws=NULL,
                                   base_w=1,
                                   ...) {
    pars <- enquos(...)
    if (xor(is.null(stanfit), is.null(PredictDiff))) {
        stop("You must specify `stanfit` or `PredictDiff` but not both.")
    }
    if (xor(is.null(stanfit), is.null(draws))) {
        stop("You must specify `stanfit` or `draws` but not both.")
    }
    if (!is.null(stanfit)) {
        PredictDiff <- GetWeightMatrixPredictor(stanfit)
        # Take the quoted variables, convert to strings, and strip the index
        # for passing to Stan directly,
        par_names <-
          lapply(dots, as_label) %>%
          sapply(function(x) { sub("\\[.*$", "", x) })
        draws <- as.matrix(stanfit, par_names)
    }
    par_means <-
        CleanGather(t(colMeans(draws))) %>%
        rename(base_mean=.value)
    join_vars <- setdiff(names(par_means), "base_mean")

    Predict <- function(w) {
        pred_mat <- PredictDiff(w=w, draws=draws, base_w=base_w)
        pred_df <-
          inner_join(
            CleanGather(pred_mat) %>%
              rename(pred_diff=.value),
            par_means,
            by=join_vars) %>%
          mutate(pred_mean=base_mean + pred_diff)
        return(pred_df)
    }
    return(Predict)
}
