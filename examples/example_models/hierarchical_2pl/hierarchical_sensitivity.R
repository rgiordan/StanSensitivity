library(rstan)
rstan_options(auto_write=TRUE)
library(ggplot2)
library(dplyr)
library(reshape2)
library(rstansensitivity)


# Set this to be the appropriate location of the repo on your computer.
example_directory <- file.path(
  Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/examples/example_models/hierarchical_2pl")


library(mirt)
sat <- key2binary(SAT12,
                  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))

sat_data <- list(I = ncol(sat),
                 J = nrow(sat),
                 N = length(sat),
                 ii = rep(1:ncol(sat), each = nrow(sat)),
                 jj = rep(1:nrow(sat), times = ncol(sat)),
                 y = as.vector(sat))

num_warmup_samples <- num_samples <- 1000
mcmc_time <- Sys.time()
model <- stan_model(file.path(example_directory, "hierarchical_2pl.stan"))
sampling_result <- sampling(model, data=sat_data, chains=1, iter=num_warmup_samples + num_samples)
mcmc_time <- Sys.time() - mcmc_time
