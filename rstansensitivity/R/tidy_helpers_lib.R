GetTidySummary <- function(stanfit, ..., spread=FALSE) {
  # Get a tidy version of a Stan summary.
  pars <- enquos(...)
  summary_mat <- t(rstan::summary(stanfit)$summary)
  tidy_summary <-
    gather_draws(summary_mat, !!!pars) %>%
    inner_join(data.frame(.draw=1:nrow(summary_mat),
                          metric=rownames(summary_mat)),
               by=".draw") %>%
    select(-.chain, -.draw, -.iteration)

  if (spread) {
    key_string <- paste(setdiff(names(tidy_summary),
                        c("metric", ".value")), collapse=" + ")
    tidy_summary <- dcast(
        tidy_summary,
        formula(sprintf("%s ~ metric", key_string)),
        value.var=".value")
  }
  return(tidy_summary)
}
