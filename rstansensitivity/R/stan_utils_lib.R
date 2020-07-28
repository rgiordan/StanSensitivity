# Return a stanfit object with no draws.  This is useful for evaluating
# log probabilities, converting between constrained and unconstrained
# parameterizations, and so on.
#' @param model A Stan model, as in the argument to \code{rstan::sampling}.
#' @param data Data, as in the argument to \code{rstan::sampling}.
#' @param init Initial parameters, as in the argument to \code{rstan::sampling}.
#' @return A \code{stanfit} object with no useful draws.
#' @export
GetDummyStanfit <- function(model, data, init=NULL) {
    if (!is.null(init)) {
        suppressWarnings(
            dummy_stanfit <-
                stan(model, data=data,
                     algorithm="Fixed_param", init=init,
                     iter=1, chains=1, refresh=0))
    } else {
        suppressWarnings(
            dummy_stanfit <-
                stan(model, data=data,
                     algorithm="Fixed_param",
                     iter=1, chains=1, refresh=0))
    }
    return(dummy_stanfit)
}


# Just a more readable shortcut for the Stan attribute.
GetParamNames <- function(model_fit) {
    model_fit@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)
}
