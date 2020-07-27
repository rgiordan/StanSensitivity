# TODO: do this better
GetGroupedLoglikMat <- function(rstan_fit, rstanarm_ij_config, df, full_loglik_mat) {
    if (rstanarm_ij_config$exchangeable_col != "") {
        stopifnot(rstanarm_ij_config$exchangeable_col %in% names(df))
        exch_col <- df[[rstanarm_ij_config$exchangeable_col]]
        loglik_mat <-
            lapply(unique(exch_col),
                   function(group_id) {
                       full_loglik_mat[, exch_col == group_id, drop=FALSE ]
                   }) %>%
            lapply(function(x) { apply(x, MARGIN=1, FUN=sum) }) %>%
            reduce(cbind)
    }
    return (loglik_mat)
}
