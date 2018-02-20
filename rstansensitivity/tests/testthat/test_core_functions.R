library(rstan)
library(rstansensitivity)
library(testthat)
rstan_options(auto_write=TRUE)

test_that("basic_functionality", {
  TypeForTest <- function(typename, varname=NULL, dim="", dist=NULL) {
  if (is.null(varname)) { 
    varname <- paste(typename, "x", sep="_")
  }
  declaration <- paste("  ", typename, dim, " ", varname, ";", sep="")
  if (is.null(dist)) {
    dist <- paste("  ", varname, " ~ normal(0, 1);", sep="")
  }
  return(data.frame(declaration=declaration, distribution=dist))
}

types_for_test <- do.call(rbind, list(
  TypeForTest("real"),
  TypeForTest("real", dim="<lower=0, upper=1>", varname="constrained_x"),
  TypeForTest("vector", dim="[2]"),
  TypeForTest("simplex", dim="[2]"),
  #TypeForTest("unit_vector", dim="[2]"), # Buggy!
  TypeForTest("ordered", dim="[2]"),
  TypeForTest("positive_ordered", dim="[2]"),
  TypeForTest("row_vector", dim="[2]"),
  TypeForTest("matrix", dim="[3,2]",
              dist="  for (i in 1:3) { for (j in 1:2) { matrix_x[i, j] ~ normal(0, 1); }}"),
  TypeForTest("corr_matrix", dim="[3]",
              dist="  corr_matrix_x ~ lkj_corr(10.0);"),
  TypeForTest("cholesky_factor_corr", dim="[3]",
              dist="  cholesky_factor_corr_x ~ lkj_corr_cholesky(10.0);"),
  TypeForTest("cov_matrix", dim="[3]",
              dist="  cov_matrix_x ~ wishart(10.0, id_mat);")
))

parameters_block <- paste(types_for_test$declaration, collapse="\n")
model_block <- paste(types_for_test$distribution, collapse="\n")

base_model_code <- paste("
data { matrix[3, 3] id_mat; }
hyperparameters {
  real theta;
}
parameters {", parameters_block, "}",
"model {", model_block, "}", sep="\n")

tmp_file <- "/tmp/test_core_functions.stan"
cat(base_model_code, file=tmp_file)
model_name <- GenerateSensitivityFromModel(tmp_file)

num_chains <- 2
model <- stan_model(GetSamplingModelFilename(model_name))
stan_data <- list(id_mat=diag(3), theta=1)
sampling_result <- sampling(model, chains=num_chains, iter=300, data=stan_data)

sens_list <- GetStanSensitivityModel(model_name, stan_data)
    
draws_array <- extract(sampling_result, permute=FALSE)
num_warmup_samples <- sampling_result@sim$warmup
num_samples <- sampling_result@sim$iter - num_warmup_samples
expect_equal(expected=num_samples, dim(draws_array)[1])
expect_equal(expected=num_chains, dim(draws_array)[2])

draws_mat <- StackChainArray(draws_array)
# Make sure there's variance so the tests are valid.
expect_true(all(apply(draws_mat, MARGIN=1, sd) > 0))

# Check the first few iterations for each chain.
for (chain in 1:num_chains) {
  for (iter in 1:3) {
    expect_equal(expected=draws_array[iter, chain, ],
                 draws_mat[(chain - 1) * num_samples + iter, ])
  }
}

# Check that get_inits is doing what it should do.
par_list_1 <- get_inits(sampling_result)[[1]]
par_list_2 <- get_inits(sampling_result)[[2]]
par_free_1 <- unconstrain_pars(sampling_result, par_list_1)
par_free_2 <- unconstrain_pars(sampling_result, par_list_2)

# Make sure they are different so the test is valid.
expect_true(sum(abs(par_free_1 - par_free_2)) > 0.01)

for (par in names(par_list_2)) {
  sens_list[[par]] <- par_list_1[[par]]
  par_list_2[[par]] <- par_list_1[[par]]  
}
par_free_2 <- unconstrain_pars(sampling_result, par_list_2)

# Now they should be equal.
expect_equal(expected=par_free_1, par_free_2, tol=1e-6)

glp_1 <- grad_log_prob(sampling_result, par_free_1)
glp_2 <- grad_log_prob(sampling_result, par_free_2)
expect_equal(glp_1, glp_2, tol=1e-6)
expect_equal(attr(glp_1, "log_prob"), attr(glp_2, "log_prob"), tol=1e-6)


# Make sure the same thing happens for the sensitivity parameters.
iter <- 20
chain <- 2
sens_par_list_1 <- sens_list$sens_par_list
par_list_1 <- get_inits(sampling_result, iter=num_warmup_samples + iter)[[chain]]
for (par in ls(par_list_1)) {
  sens_par_list_1[[par]] <- par_list_1[[par]]
}
par_free_1 <- unconstrain_pars(sampling_result, par_list_1)
sens_par_free_1 <- unconstrain_pars(sens_list$model_sens_fit, sens_par_list_1)
glp_1 <- grad_log_prob(sampling_result, par_free_1)
glp_2 <- grad_log_prob(sens_list$model_sens_fit, sens_par_free_1)
lp_col <- dim(draws_array)[3]
expect_equal(glp_1[1:length(glp_1)], glp_2[1:length(glp_1)], tol=1e-6)
expect_equal(attr(glp_1, "log_prob"), attr(glp_2, "log_prob"), tol=1e-6)
expect_equal(attr(glp_1, "log_prob"), draws_array[iter, chain, lp_col], tol=1e-6)


##################
# Check EvaluateAtDraws.

sens_par_list <- get_inits(sens_list$model_sens_params)[[1]]
sens_par_list$theta <- stan_data$theta

draws_list <- EvaluateAtDraws(
  sampling_result, sens_list, sens_par_list, compute_grads=TRUE)
lp_col <- dim(draws_array)[3]
expect_equal(draws_list$lp_vec, as.numeric(draws_array[,, lp_col]), tol=1e-6)

})

