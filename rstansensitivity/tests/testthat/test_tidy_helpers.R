#!/usr/bin/env Rscript

library(devtools)
library(rstansensitivity)
library(testthat)
rstan_options(auto_write=TRUE)

devtools::load_all()
context("rstansensitivity")


test_that("TidyColumnWorks", {
    parnames <- c(sprintf("a[%d]", 1:3), "b", sprintf("c[%d]", 1:5), "d")

    df <- data.frame(par=parnames, val=runif(length(parnames)))
    col <- "par"
    pars <- list(quote(a[i]))

    tidy_df <- TidyColumn(df, "par", a[d], variable_name="hyperparameter")

    for (this_d in 1:3) {
        df_row <- filter(df, par == sprintf("a[%d]", this_d))
        tidy_df_row <- filter(tidy_df, hyperparameter == "a", d == this_d)
        expect_equal(df_row$val, tidy_df_row$val, tol=1e-12)
        expect_equal(as.character(df_row$par), as.character(tidy_df_row$par))
    }
})


test_that("GatherRowNamedMatrixWorks", {
    pars <- c("a", "b", "cde")
    mat <- matrix(runif(12), nrow=3, ncol=4,
                  dimnames=list(pars, c("x[1]", "x[2]", "y", "z")))
    df <- GatherRowNamedMatrix(mat, x[i], y)
    standard_cols <- c("i", ".variable", ".value")
    expect_equal(gather_draws(mat, x[i], y)[standard_cols], df[standard_cols])
    expect_equal(pars, unique(as.character(df$.rowname)))
})
