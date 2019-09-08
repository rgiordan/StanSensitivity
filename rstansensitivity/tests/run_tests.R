#!/usr/bin/env Rscript
library(testthat)

cat("\nNote: the first time the tests are run will be slow due to the need",
    "to compile the Stan models.  Subsequent runs will be much faster.\n\n")

# I'm not 100% sure why or whether this separate file is necessary.
test_dir(".")
